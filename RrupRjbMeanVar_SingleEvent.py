#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import sys
import os
import os.path
import datetime
try:
    import configparser
except ImportError:
    import ConfigParser as configparser 
import time as time
import numpy as np
from scipy.stats import norm

#####################################################
# This function loops over M and R, calling EpsilonDipDxDyTheta
# and writes out progress information.
#####################################################
def RrupMeanVar_MR(ndip=19, mindip=0.0, maxdip=np.pi/2.0,
                   min_seis_depth=0, max_seis_depth=20,
                   rup_dim_model='WC94', mech="A", AR=1,
                   ntheta=73, nxny=100, zhyp=0, bytheta=False,
                   neps=10, trunc=3,
                   NP=1, iP=0, rjb_filename='junk', rrup_filename='junk',
                   M=6.0, Repi=100):

    nepi = np.size(Repi)

    if NP != 1:
        ii = range(iP, nepi, NP)
        Repi = Repi[ii]
        nepi = np.size(Repi)

    if bytheta is True:
        Rrup_ratio = np.zeros((nepi, ntheta))
        Rrup_variance = np.zeros((nepi, ntheta))
        Rjb_ratio = np.zeros((nepi, ntheta))
        Rjb_variance = np.zeros((nepi, ntheta))
    else:
        Rrup_ratio = np.zeros((nepi, 1))
        Rrup_variance = np.zeros((nepi, 1))
        Rjb_ratio = np.zeros((nepi, 1))
        Rjb_variance = np.zeros((nepi, 1))

    for i in range(0, nepi):
        if bytheta is True:
            theta = np.linspace(0, np.pi*2, ntheta)  ## in rad
            for j in range(0, ntheta):
                rrup_var, rrup_avg, rjb_var, rjb_avg = \
                    EpsilonDipDxDyTheta(M, Repi[i],
                                        ndip=ndip, mindip=mindip, maxdip=maxdip,
                                        min_seis_depth=min_seis_depth,
                                        max_seis_depth=max_seis_depth,
                                        rup_dim_model=rup_dim_model, mech=mech, AR=AR,
                                        theta=np.array([theta[j]]), bytheta=bytheta,
                                        ntheta=1, nxny=nxny, zhyp=zhyp,
                                        neps=neps, trunc=trunc)
                Rrup_ratio[i, j] = rrup_avg / Repi[i]
                Rrup_variance[i, j] = rrup_var
                Rjb_ratio[i, j] = rjb_avg / Repi[i]
                Rjb_variance[i, j] = rjb_var
                print("        Proc %d j=%d of %d %s" % (iP, j+1, ntheta,
                    datetime.datetime.now().isoformat()))
                print("Proc %d done with %d of %d distances %s" % (iP, i+1, nepi,
                    datetime.datetime.now().isoformat()))
        else:
            rrup_var, rrup_avg, rjb_var, rjb_avg = \
                EpsilonDipDxDyTheta(M, Repi[i],
                                    ndip=ndip, mindip=mindip, maxdip=maxdip,
                                    min_seis_depth=min_seis_depth,
                                    max_seis_depth=max_seis_depth,
                                    rup_dim_model=rup_dim_model, mech=mech, AR=AR,
                                    theta=0, bytheta=bytheta,
                                    ntheta=ntheta, nxny=nxny, zhyp=zhyp,
                                    neps=neps, trunc=trunc)
            Rrup_ratio[i, 0] = rrup_avg / Repi[i]
            Rrup_variance[i, 0] = rrup_var
            Rjb_ratio[i, 0] = rjb_avg / Repi[i]
            Rjb_variance[i, 0] = rjb_var
            print("        Proc %d j=%d of %d %s" % (iP, 1, 1,
                datetime.datetime.now().isoformat()))
            print("Proc %d done with %d of %d distances %s" % (iP, i+1, nepi,
                datetime.datetime.now().isoformat()))
    fr_rrup = open('%sRatios_%02d.csv' % (rrup_filename, iP), 'w')
    fv_rrup = open('%sVar_%02d.csv' % (rrup_filename, iP), 'w')
    fr_rjb = open('%sRatios_%02d.csv' % (rjb_filename, iP), 'w')
    fv_rjb = open('%sVar_%02d.csv' % (rjb_filename, iP), 'w')

    for i in range(0, nepi):
        fr_rrup.write("%f," % Repi[i])
        fv_rrup.write("%f," % Repi[i])
        fr_rjb.write("%f," % Repi[i])
        fv_rjb.write("%f," % Repi[i])
        if bytheta is True:
            for j in range(0, ntheta):
                fr_rrup.write("%f" % Rrup_ratio[i,j])
                fv_rrup.write("%f" % Rrup_variance[i,j])
                fr_rjb.write("%f" % Rjb_ratio[i,j])
                fv_rjb.write("%f" % Rjb_variance[i,j])
                if j < ntheta - 1:
                    fr_rrup.write(",")
                    fv_rrup.write(",")
                    fr_rjb.write(",")
                    fv_rjb.write(",")
            fr_rrup.write("\n")
            fv_rrup.write("\n")
            fr_rjb.write("\n")
            fv_rjb.write("\n")
        else:
            fr_rrup.write("%f\n" % Rrup_ratio[i,0])
            fv_rrup.write("%f\n" % Rrup_variance[i,0])
            fr_rjb.write("%f\n" % Rjb_ratio[i,0])
            fv_rjb.write("%f\n" % Rjb_variance[i,0])

    fr_rrup.close()
    fv_rrup.close()
    fr_rjb.close()
    fv_rjb.close()
    sys.exit(0)

#####################################################
# This function evaluates the Rrup mean and var integral
# for a single M/R pair, looping over:
#    - dip
#    - dx, dy (location of hypocenter on fault)
#    - theta (angle to fault)
#    - epsilon (dummy variable for L/W/A integration)
#####################################################
def EpsilonDipDxDyTheta(M, Repi, ndip=19, mindip=0.0, maxdip=np.pi/2.0,
                        min_seis_depth=0, max_seis_depth=20,
                        rup_dim_model='WC94', mech="A", AR=1,
                        theta=0, bytheta=False, ntheta=73, nxny=50, zhyp=0,
                        neps=10, trunc=3):
    ## Same as before, but position of P is random so have to integrate from
    ## theta = 0 to 2pi and then for dx = 0 to L and dy = 0 to SW
    ## and dip from 0 to pi/2
    if ndip != 1:
        dip = np.linspace(mindip, maxdip, ndip)
        ddip = dip[1] - dip[0]
        dipnorm = 1.0 / (maxdip - mindip)
    else:
        dip = np.array([mindip])

    ## Add integration across length and width.
    ## Need to assume a truncation level for normal distribution
    eps = np.linspace(-trunc, trunc, neps+1)
    epsmid = 0.5*(eps[1:] + eps[:-1])
    peps = (norm.cdf(eps[1:]) - norm.cdf(eps[:-1]))

    #
    # define delta epsilons to normalize probabilities
    #
    d_eps = 2 * trunc / neps
    epsfac = np.trapz(peps, dx=d_eps)

    if neps > 1:
        peps /= epsfac

    if rup_dim_model == 'WC94':
        ## Use mech to get either M-A or (M-W) and (M-R) from Wells and Coppersmith.
        ## Note: these relations are not specific to active environments, but there
        ##       are specific equations for stable that we can use in that case. But
        ##       it wouldn't be 'wrong' to use these for stable.
        if mech == "SS":
            siga = 0.22
            A = np.power(10, -3.42 + 0.90*M + siga * epsmid)
        if mech == "R":
            siga = 0.26
            A = np.power(10, -3.99 + 0.98*M + siga * epsmid)
        if mech == "N":
            siga = 0.22
            A = np.power(10, -2.78 + 0.82*M + siga * epsmid)
        if mech == "A":
            siga = 0.24
            A = np.power(10, -3.49 + 0.91*M + siga * epsmid)
    else:
        # Somerville (2014), used in CEUS-SSC
        siga = 0.3
        A = np.power(10, M - 4.25 + siga * epsmid)

    if bytheta is False:
        theta = np.linspace(0, 2*np.pi, ntheta)  ## in rad
        dt = theta[1] - theta[0]


    ## origin defined at o:
    ##
    ##  + +------+
    ##  | |      |
    ##  L |      |
    ##  | |      |
    ##  + o------+
    ##    +--SW--+
    ##
    one_over_2pi = 1.0/2.0/np.pi

    RrupIntegrand_a = np.zeros(neps) + np.nan
    RrupIntegrand_a2 = np.zeros(neps) + np.nan
    RrupIntegrand_d = np.zeros(ndip)
    RrupIntegrand_d2 = np.zeros(ndip)
    RrupIntegrand_y = np.zeros(nxny)
    RrupIntegrand_y2 = np.zeros(nxny)
    RrupIntegrand_x = np.zeros(nxny)
    RrupIntegrand_x2 = np.zeros(nxny)

    RjbIntegrand_a = np.zeros(neps) + np.nan
    RjbIntegrand_a2 = np.zeros(neps) + np.nan
    RjbIntegrand_d = np.zeros(ndip)
    RjbIntegrand_d2 = np.zeros(ndip)
    RjbIntegrand_y = np.zeros(nxny)
    RjbIntegrand_y2 = np.zeros(nxny)
    RjbIntegrand_x = np.zeros(nxny)
    RjbIntegrand_x2 = np.zeros(nxny)

    Rjb = np.zeros(ntheta)
    Rrup = np.zeros(ntheta)
    Rrupp = np.zeros(ntheta)
    Ry = np.zeros(ntheta)
    nx = nxny
    ny = nxny


    for m in range(0, neps): # area
        W = np.sqrt(A[m]/AR)
        for k in range(0, ndip): # dip
            if np.allclose(dip[k], 0) is False:
                ZW = W * np.sin(dip[k]) # vertical projection of W
                # Overwrite W if it extends too far and recompute SW, x, dx
                Ztor = np.max(np.array([zhyp - ZW, min_seis_depth]))
                Zbor = np.min(np.array([zhyp + ZW, max_seis_depth]))
                W = (Zbor - Ztor)/np.sin(dip[k])
            else:
                Ztor = zhyp

            L = A[m]/W
            SW = W * np.cos(dip[k])
            x = np.linspace(0, SW, nx)
            dx = x[1] - x[0]
            one_over_L = 1.0 / L
            one_over_sw = 1.0 / SW
            y = np.linspace(0, L, ny)
            dy = y[1] - y[0]
            for i in range(0, nx): # x
                xj = x[i] + Repi * np.cos(theta)

                xltz = xj < 0
                xgez_and_xlesw = (xj >= 0) & (xj <= SW)
                xgtsw = xj > SW
                c1x = xltz
                c2x = xgez_and_xlesw
                c3x = xgtsw
                c4x = xltz
                c5x = xgez_and_xlesw
                c6x = xgtsw
                c7x = xltz
                c8x = xgez_and_xlesw
                c9x = xgtsw
                for j in range(0, ny):
                    yi = y[j] + Repi * np.sin(theta)

                    cca = yi > L
                    ccb = (yi >= 0) & (yi <= L)
                    ccc = yi < 0
                    c1 = c1x & cca
                    c2 = c2x & cca
                    c3 = c3x & cca
                    c4 = c4x & ccb
                    c5 = c5x & ccb
                    c6 = c6x & ccb
                    c7 = c7x & ccc
                    c8 = c8x & ccc
                    c9 = c9x & ccc
                    xx = xj[c1]
                    yy = yi[c1] - L
                    Rjb[c1] = np.sqrt(xx*xx + yy*yy)
                    Rjb[c2] = yi[c2] - L
                    xx = xj[c3] - SW
                    yy = yi[c3] - L
                    Rjb[c3] = np.sqrt(xx*xx + yy*yy)
                    Rjb[c4] = np.abs(xj[c4])
                    Rjb[c5] = 0
                    Rjb[c6] = xj[c6] - SW
                    xx = xj[c7]
                    yy = yi[c7]
                    Rjb[c7] = np.sqrt(xx*xx + yy*yy)
                    Rjb[c8] = np.abs(yi[c8])
                    xx = xj[c9] - SW
                    yy = yi[c9]
                    Rjb[c9] = np.sqrt(xx*xx + yy*yy)
                    #################################
                    # Compute Rx and Ry
                    Rx = xj
                    Ry[c1 | c2 | c3] = yi[c1 | c2 | c3] - L
                    Ry[c4 | c5 | c6] = 0
                    Ry[c7 | c8 | c9] = np.abs(yi[c7 | c8 | c9])
                    #################################
                    # Compute Rrup prime, then Rrup
                    # using Kaklamanos eqns
                    if dip[k] == np.pi/2:
                        Rrup = np.sqrt(Rjb**2 + Ztor**2)
                    else:
                        tmp = (Ztor * np.tan(dip[k]))
                        r1 = Rx < tmp
                        Rrupp[r1] = np.sqrt(Rx[r1]**2 + Ztor**2)
                        r2 = (tmp <= Rx) & (Rx <= (tmp + W/np.cos(dip[k])))
                        Rrupp[r2] = Rx[r2]*np.sin(dip[k]) + Ztor*np.cos(dip[k])
                        r3 = Rx > (tmp + W/np.cos(dip[k]))
                        Rrupp[r3] = np.sqrt((Rx[r3] - W*np.cos(dip[k]))**2 + \
                                    (Ztor + W*np.sin(dip[k]))**2)

                    Rrup = np.sqrt(Rrupp**2 + Ry**2)
                    Rrup2 = Rrup * Rrup
                    Rjb2 = Rjb * Rjb
                    if bytheta is False:
                        RrupIntegrand_y[j] = one_over_2pi * np.trapz(Rrup, dx=dt)
                        RrupIntegrand_y2[j] = one_over_2pi * np.trapz(Rrup2, dx=dt)
                        RjbIntegrand_y[j] = one_over_2pi * np.trapz(Rjb, dx=dt)
                        RjbIntegrand_y2[j] = one_over_2pi * np.trapz(Rjb2, dx=dt)
                    else:
                        RrupIntegrand_y[j] = Rrup
                        RrupIntegrand_y2[j] = Rrup2
                        RjbIntegrand_y[j] = Rjb
                        RjbIntegrand_y2[j] = Rjb2

                RrupIntegrand_x[i] = one_over_L * np.trapz(RrupIntegrand_y, dx=dy)
                RrupIntegrand_x2[i] = one_over_L * np.trapz(RrupIntegrand_y2, dx=dy)
                RjbIntegrand_x[i] = one_over_L * np.trapz(RjbIntegrand_y, dx=dy)
                RjbIntegrand_x2[i] = one_over_L * np.trapz(RjbIntegrand_y2, dx=dy)

            RrupIntegrand_d[k] = one_over_sw * np.trapz(RrupIntegrand_x, dx=dx)
            RrupIntegrand_d2[k] = one_over_sw * np.trapz(RrupIntegrand_x2, dx=dx)
            RjbIntegrand_d[k] = one_over_sw * np.trapz(RjbIntegrand_x, dx=dx)
            RjbIntegrand_d2[k] = one_over_sw * np.trapz(RjbIntegrand_x2, dx=dx)
        if ndip == 1:
            RrupIntegrand_a[m] = RrupIntegrand_d[0]
            RrupIntegrand_a2[m] = RrupIntegrand_d2[0]
            RjbIntegrand_a[m] = RjbIntegrand_d[0]
            RjbIntegrand_a2[m] = RjbIntegrand_d2[0]
        else:
            RrupIntegrand_a[m] = dipnorm * np.trapz(RrupIntegrand_d, dx=ddip)
            RrupIntegrand_a2[m] = dipnorm * np.trapz(RrupIntegrand_d2, dx=ddip)
            RjbIntegrand_a[m] = dipnorm * np.trapz(RjbIntegrand_d, dx=ddip)
            RjbIntegrand_a2[m] = dipnorm * np.trapz(RjbIntegrand_d2, dx=ddip)
    if neps == 1:
        Rrup_avg = RrupIntegrand_a[0]
        Rrup_var = RrupIntegrand_a2[0] - Rrup_avg**2
        Rjb_avg = RjbIntegrand_a[0]
        Rjb_var = RjbIntegrand_a2[0] - Rjb_avg**2
    else:
        Rrup_avg = np.trapz(peps*RrupIntegrand_a, dx=d_eps)
        Rrup_var = np.trapz(peps*RrupIntegrand_a2, dx=d_eps) - Rrup_avg**2
        Rjb_avg = np.trapz(peps*RjbIntegrand_a, dx=d_eps)
        Rjb_var = np.trapz(peps*RjbIntegrand_a2, dx=d_eps) - Rjb_avg**2

    return Rrup_var, Rrup_avg, Rjb_var, Rjb_avg

def main():

    start_datetime = datetime.datetime.now().isoformat()

    config_file = sys.argv[-1]
    if not os.path.exists(config_file):
        raise Exception("Config file %s doesn't exist" % (config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)
    section = Config.sections()[0]

    NP = Config.getint(section, 'NP')
    datadir = Config.get(section, 'datadir')
    M = Config.getfloat(section, 'M')
    zhyp = Config.getfloat(section, 'zhyp')
    bytheta = Config.getboolean(section, 'bytheta')
    rup_dim_model = Config.get(section, 'rup_dim_model')
    mech = Config.get(section, 'mech')
    AR = Config.getfloat(section, 'AR')
    ndip = Config.getint(section, 'ndip')
    mindip_deg = Config.getfloat(section, 'mindip_deg')
    maxdip_deg = Config.getfloat(section, 'maxdip_deg')
    ntheta = Config.getint(section, 'ntheta')
    nxny = Config.getint(section, 'nxny')
    neps = Config.getint(section, 'neps')
    trunc = Config.getfloat(section, 'trunc')
    minepi = Config.getfloat(section, 'minepi')
    maxepi = Config.getfloat(section, 'maxepi')
    nepi = Config.getint(section, 'nepi')
    min_seis_depth = Config.getfloat(section, 'min_seis_depth')
    max_seis_depth = Config.getfloat(section, 'max_seis_depth')

    if os.path.isdir(datadir) is False:
        os.mkdir(datadir)

    theta_string = ''
    if bytheta is True:
        theta_string = '_bytheta'

    rjb_filename = os.path.join(datadir, "Rjb%s" % (theta_string))
    rrup_filename = os.path.join(datadir, "Rrup%s" % (theta_string))

    mindip = mindip_deg * np.pi / 180.0
    maxdip = maxdip_deg * np.pi / 180.0
    Repi = np.logspace(np.log10(minepi), np.log10(maxepi), nepi)
    #Repi = np.hstack((np.logspace(np.log10(0.1), np.log10(1), 6)[0:-1],
    #                  np.logspace(np.log10(1), np.log10(1000), 31)))

    for iP in range(0, NP):
        if os.fork() == 0:
            RrupMeanVar_MR(rup_dim_model=rup_dim_model, mech=mech, AR=AR,
                           ndip=ndip, mindip=mindip, maxdip=maxdip,
                           min_seis_depth=min_seis_depth, max_seis_depth=max_seis_depth,
                           ntheta=ntheta, nxny=nxny, zhyp=zhyp, bytheta=bytheta, neps=neps,
                           trunc=trunc, NP=NP, iP=iP, rjb_filename=rjb_filename,
                           rrup_filename=rrup_filename, M=M, Repi=Repi)

    for iP in range(0, NP):
        pid, status = os.waitpid(-1, 0)

    fd = [None] * 4
    fd[0] = open('%s_Ratios.csv' % (rjb_filename), 'w')
    fd[1] = open('%s_Var.csv' % (rjb_filename), 'w')
    fd[2] = open('%s_Ratios.csv' % (rrup_filename), 'w')
    fd[3] = open('%s_Var.csv' % (rrup_filename), 'w')

    for i in range(0, 4):
        fd[i].write('# Program: %s\n' % sys.argv[0])
        fd[i].write('# Config file: %s\n' % config_file)
        fd[i].write('# Process start: %s\n' % start_datetime)
        fd[i].write('# Process finish: %s\n' % datetime.datetime.now().isoformat())
        fd[i].write('# M = %f, zhyp = %f, bytheta = %s, rup_dim_model = %s, mech = %s, '
                    'AR = %s, ndip = %d, mindip = %f, maxdip = %f\n' \
                    % (M, zhyp, bytheta, rup_dim_model, mech, AR, ndip, mindip, maxdip))
        fd[i].write('# ntheta = %d, nxny = %d, neps = %d, trunc = %f, min_seis_depth = %f, '
                    'max_seis_depth = %f\n'
                    % (ntheta, nxny, neps, trunc, min_seis_depth, max_seis_depth))
        fd[i].write('"Repi_km",')

    if bytheta is True:
        theta = np.linspace(0, 2 * np.pi, ntheta)
        for i in range(0, 4):
            for j in range(0, ntheta):
                fd[i].write('"%g"' % (theta[j]))
                if j < ntheta - 1:
                    fd[i].write(',')
            fd[i].write("\n")
    else:
        fd[0].write('"R%g"\n' % M)
        fd[1].write('"V%g"\n' % M)
        fd[2].write('"R%g"\n' % M)
        fd[3].write('"V%g"\n' % M)

    rjb_rat_file = [None] * NP
    rjb_var_file = [None] * NP
    rrup_rat_file = [None] * NP
    rrup_var_file = [None] * NP
    rjb_frs = [None] * NP
    rjb_fvs = [None] * NP
    rrup_frs = [None] * NP
    rrup_fvs = [None] * NP
    for iP in range(0, NP):
        rjb_rat_file[iP] = '%sRatios_%02d.csv' % (rjb_filename, iP)
        rjb_var_file[iP] = '%sVar_%02d.csv' % (rjb_filename, iP)
        rrup_rat_file[iP] = '%sRatios_%02d.csv' % (rrup_filename, iP)
        rrup_var_file[iP] = '%sVar_%02d.csv' % (rrup_filename, iP)
        rjb_frs[iP] = open(rjb_rat_file[iP], 'r')
        rjb_fvs[iP] = open(rjb_var_file[iP], 'r')
        rrup_frs[iP] = open(rrup_rat_file[iP], 'r')
        rrup_fvs[iP] = open(rrup_var_file[iP], 'r')

    for line in rjb_frs[0]:
        fd[0].write(line)
        line = rjb_fvs[0].readline()
        fd[1].write(line)
        line = rrup_frs[0].readline()
        fd[2].write(line)
        line = rrup_fvs[0].readline()
        fd[3].write(line)
        for iP in range(1, NP):
            line = rjb_frs[iP].readline()
            if line:
                fd[0].write(line)
            line = rjb_fvs[iP].readline()
            if line:
                fd[1].write(line)
            line = rrup_frs[iP].readline()
            if line:
                fd[2].write(line)
            line = rrup_fvs[iP].readline()
            if line:
                fd[3].write(line)

    for iP in range(0, NP):
        rjb_frs[iP].close()
        rjb_fvs[iP].close()
        rrup_frs[iP].close()
        rrup_fvs[iP].close()
        os.unlink(rjb_rat_file[iP])
        os.unlink(rjb_var_file[iP])
        os.unlink(rrup_rat_file[iP])
        os.unlink(rrup_var_file[iP])

    for i in range(0, 4):
        fd[i].close()

if __name__ == "__main__":

    start_time = time.time()

    main()

    print("Total execution time %f seconds." % (time.time() - start_time))
    sys.exit(0)

