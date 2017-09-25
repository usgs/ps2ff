#!/usr/bin/env python
"""
This is a test.
"""

#
# NOTICE: This is experimental and incomplete code that is not part 
# of the reviewed and documented product found in this repository. 
# It is currently under development and is not fit to be used for 
# any purpose.
#

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

from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib.imt import PGA, PGV, SA
cy14 = ChiouYoungs2014()
coef = cy14.COEFFS.non_sa_coeffs[PGA()]

#####################################################
# This function loops over M and R, calling EpsilonDipDxDyTheta
# and writes out progress information.
#####################################################
def MeanVar_MR(ndip=19, mindip=0.0, maxdip=np.pi/2.0,
               min_seis_depth=0, max_seis_depth=20,
               rup_dim_model='WC94', mech="A", LW=True, AR=1,
               ntheta=73, nxny=100, nz=20,
               neps=10, trunc=3,
               NP=1, iP=0, filename='Cy14Hw-',
               M=6.0, Repi=100):

    nepi = np.size(Repi)

    if NP != 1:
        ii = range(iP, nepi, NP)
        Repi = Repi[ii]
        nepi = np.size(Repi)

    nmag = np.size(M)
    ratio = np.zeros((nepi, nmag))
    variance = np.zeros((nepi, nmag))
    for i in range(0, nepi):
        for j in range(0, nmag):
            var, avg = EpsilonDipDxDyTheta(M[j], Repi[i],
                                           ndip=ndip, mindip=mindip, maxdip=maxdip,
                                           min_seis_depth=min_seis_depth,
                                           max_seis_depth=max_seis_depth,
                                           rup_dim_model=rup_dim_model, mech=mech,
                                           LW=LW, AR=AR,
                                           ntheta=ntheta, nxny=nxny, nz=nz,
                                           neps=neps, trunc=trunc)
            ratio[i,j] = avg
            variance[i,j] = var
            print("        Proc %d j=%d of %d %s" % (iP, j+1, nmag,
                datetime.datetime.now().isoformat()))
        print("Proc %d done with %d of %d distances %s" % (iP, i+1, nepi,
            datetime.datetime.now().isoformat()))

    fr = open('%sHw_%02d.csv' % (filename, iP), 'w')
    fv = open('%sHwVar_%02d.csv' % (filename, iP), 'w')

    for i in range(0, nepi):
        fr.write("%f," % Repi[i])
        fv.write("%f," % Repi[i])
        for j in range(0, nmag):
            fr.write("%f" % ratio[i,j])
            fv.write("%f" % variance[i,j])
            if j < nmag - 1:
                fr.write(",")
                fv.write(",")
        fr.write("\n")
        fv.write("\n")

    fr.close()
    fv.close()
    sys.exit(0)

#####################################################
# This function evaluates the mean and var integral
# for a single M/R pair, looping over:
#    - dip
#    - dx, dy (location of hypocenter on fault)
#    - theta (angle to fault)
#    - epsilon (dummy variable for L/W/A integration)
#####################################################
def EpsilonDipDxDyTheta(M, Repi, ndip=19, mindip=0.0, maxdip=np.pi/2.0,
                        min_seis_depth=0,  max_seis_depth=20,
                        rup_dim_model='WC94', mech="A", LW=True, AR=1,
                        ntheta=73, nxny=50, nz=20,
                        neps=10, trunc=3):
    max_seis_thickness = max_seis_depth - min_seis_depth
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
            if LW is True:
                sigl = 0.15
                L = 10**(-2.57 + 0.62 * M + sigl * epsmid)
                sigw = 0.14
                W = 10**(-0.76 + 0.27 * M + sigw * epsmid)
                nl = neps
                nw = neps
            else:
                siga = 0.22
                A = np.power(10, -3.42 + 0.90*M + siga * epsmid.reshape(1,-1))
        if mech == "R":
            if LW is True:
                sigl = 0.16
                L = 10**(-2.42 + 0.58 * M + sigl * epsmid)
                sigw = 0.15
                W = 10**(-1.61 + 0.41 * M + sigw * epsmid)
                nl = neps
                nw = neps
            else:
                siga = 0.26
                A = np.power(10, -3.99 + 0.98*M + siga * epsmid.reshape(1,-1))
        if mech == "N":
            if LW is True:
                sigl = 0.17
                L = 10**(-1.88 + 0.50 * M) + sigl * epsmid
                sigw = 0.12
                W = 10**(-1.14 + 0.35 * M + sigw * epsmid)
                nl = neps
                nw = neps
            else:
                siga = 0.22
                A = np.power(10, -2.78 + 0.82*M + siga * epsmid.reshape(1,-1))
        if mech == "A":
            if LW is True:
                sigl = 0.16
                L = 10**(-2.44 + 0.59 * M + sigl * epsmid)
                sigw = 0.15
                W = 10**(-1.01 + 0.32 * M + sigw * epsmid)
                nl = neps
                nw = neps
            else:
                siga = 0.24
                A = np.power(10, -3.49 + 0.91*M + siga * epsmid.reshape(1,-1))
    else:
        siga = 0.3
        A = np.power(10, M - 4.25 + siga * epsmid.reshape(1,-1))

    if rup_dim_model != 'WC94' or LW is False:
        # Trick the L and W loops to handle a one to one mapping
        # between L and W, and each L/W pair corresponds to a
        # single M looping over epsilon; specify length and
        # width constrained by seismogenic depth

        W_mat = np.tile(np.sqrt(A/AR), (ndip, 1))
        sindip = np.tile(np.sin(dip).reshape(-1,1), (1,neps))
        rup_z = W_mat * sindip
        indxx = rup_z > max_seis_thickness
        W_mat[indxx] = max_seis_thickness / sindip[indxx]
        L_mat = np.tile(A, (ndip, 1)) / W_mat
        nl = 1
        nw = neps

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

    integrand = type('integrand', (object,), {})()
    integrand.w = np.zeros(nw) + np.nan
    integrand.w2 = np.zeros(nw) + np.nan
    integrand.l = np.zeros(nl)
    integrand.l2 = np.zeros(nl)
    integrand.d = np.zeros(ndip)
    integrand.d2 = np.zeros(ndip)
    integrand.z = np.zeros(nz)
    integrand.z2 = np.zeros(nz)
    Rjb = np.zeros(ntheta)
    Rrup = np.zeros(ntheta)
    Rrupp = np.zeros(ntheta)
    Ry = np.zeros(ntheta)

#    mindx = Repi / M / 2.0

    # This is the loop over W
    for w in range(0, nw):
        # This is the loop over L
        for l in range(0, nl):
            for k in range(0, ndip):
                if rup_dim_model == 'WC94' and LW is True:
                    ll = l
                else:
                    W = W_mat[k,:]
                    L = L_mat[k,:]
                    ll = w
                one_over_ll = 1.0 / L[ll]

                # Since we still use linspace, dy won't be exactly minx.
                # Also put a hard bound on minimum ny to be 2 so that dy makes sense.
                ny = nxny
                y = np.linspace(0, L[ll], ny)
                dy = y[1] - y[0]

                integrand.y = np.zeros(ny)
                integrand.y2 = np.zeros(ny)
                ## Fault width projected to surface:
                if np.allclose(np.cos(dip[k]), 0):
                    SW = 0
                    nx = 1
                    x = np.linspace(0, SW, nx)
                    dx = 0
                else:
                    SW = W[w] * np.cos(dip[k])
                    one_over_sw = 1.0 / SW

                    # Since we still use linspace, dx won't be exactly minx.
                    # Also put a hard bound on minimum nx to be 2 so that dx makes sense.
                    nx = nxny
                    x = np.linspace(0, SW, nx)
                    dx = x[1] - x[0]
                # Calclate range of Ztor
                ZtorMax = max_seis_depth - W[w]*np.sin(dip[k])
                if np.allclose(ZtorMax, 0) or nz == 1:
                    nz = 1
                    dz = 0
                    Ztor = np.linspace(0, ZtorMax, nz)
                else:
                    one_over_ztormax = 1.0/ZtorMax
                    Ztor = np.linspace(0, ZtorMax, nz)
                    dz = Ztor[1]-Ztor[0]


                integrand.x = np.zeros(nx)
                integrand.x2 = np.zeros(nx)
                for z in range(0, nz):
                    for i in range(0, nx):
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

                            cca = yi > L[ll]
                            ccb = (yi >= 0) & (yi <= L[ll])
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
                            # The original version. The above two chunks should
                            # implement this
                            # c1 = (xj <  0) &              (yi > L[l])
                            # c2 = (xj >= 0) & (xj <= SW) & (yi > L[l])
                            # c3 = (xj > SW) &              (yi > L[l])
                            # c4 = (xj <  0) &              (yi >= 0) & (yi <= L[l])
                            # c5 = (xj >= 0) & (xj <= SW) & (yi >= 0) & (yi <= L[l])
                            # c6 = (xj > SW) &              (yi >= 0) & (yi <= L[l])
                            # c7 = (xj <  0) &              (yi < 0)
                            # c8 = (xj >= 0) & (xj <= SW) & (yi < 0)
                            # c9 = (xj > SW) &              (yi < 0)

                            xx = xj[c1]
                            yy = yi[c1] - L[ll]
                            Rjb[c1] = np.sqrt(xx*xx + yy*yy)
                            Rjb[c2] = yi[c2] - L[ll]
                            xx = xj[c3] - SW
                            yy = yi[c3] - L[ll]
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
                            Ry[c1 | c2 | c3] = yi[c1 | c2 | c3] - L[ll]
                            Ry[c4 | c5 | c6] = 0
                            Ry[c7 | c8 | c9] = np.abs(yi[c7 | c8 | c9])
                            #################################
                            # Compute Rrup prime, then Rrup
                            # using Kaklamanos eqns
                            if dip[k] == np.pi/2:
                                Rrup = np.sqrt(Rjb**2 + Ztor[z]**2)
                            else:
                                tmp = (Ztor[z] * np.tan(dip[k]))
                                r1 = Rx < tmp
                                Rrupp[r1] = np.sqrt(Rx[r1]**2 + Ztor[z]**2)
                                r2 = (tmp <= Rx) & (Rx <= (tmp + W[w]/np.cos(dip[k])))
                                Rrupp[r2] = Rx[r2]*np.sin(dip[k]) + Ztor[z]*np.cos(dip[k])
                                r3 = Rx > (tmp + W[w]/np.cos(dip[k]))
                                Rrupp[r3] = np.sqrt((Rx[r3] - W[w]*np.cos(dip[k]))**2 + (Ztor[z] + W[w]*np.sin(dip[k]))**2)

                            Rrup = np.sqrt(Rrupp**2 + Ry**2)
                            Rrup2 = Rrup * Rrup

                            Fhw = np.array(Rx >= 0, dtype=float)
                            hw1 = coef['c9'] * Fhw * np.cos(dip[k])
                            hw2 = coef['c9a'] + (1 - coef['c9a']) * np.tanh(Rx/coef['c9b'])
                            hw3 = 1 - np.sqrt(Rjb**2 + Ztor[z]**2)/(Rrup + 1)
                            hw = hw1 * hw2 * hw3
                            hw2 = hw * hw

                            integrand.y[j] = one_over_2pi * np.trapz(hw, dx=dt)
                            integrand.y2[j] = one_over_2pi * np.trapz(hw2, dx=dt)

                        integrand.x[i] = one_over_ll * np.trapz(integrand.y, dx=dy)
                        integrand.x2[i] = one_over_ll * np.trapz(integrand.y2, dx=dy)

                    if dx == 0:
                        integrand.z[z] = integrand.x[0]
                        integrand.z2[z] = integrand.x2[0]
                    else:
                        integrand.z[z] = one_over_sw * np.trapz(integrand.x, dx=dx)
                        integrand.z2[z] = one_over_sw * np.trapz(integrand.x2, dx=dx)
                # end z loop
                if dz == 0:
                    integrand.d[k] = integrand.z[0]
                    integrand.d2[k] = integrand.z2[0]
                else:
                    integrand.d[k] = one_over_ztormax * np.trapz(integrand.z, dx=dz)
                    integrand.d2[k] = one_over_ztormax * np.trapz(integrand.z2, dx=dz)
            if ndip == 1:
                integrand.l[l] = integrand.d[0]
                integrand.l2[l] = integrand.d2[0]
            else:
                integrand.l[l] = dipnorm * np.trapz(integrand.d, dx=ddip)
                integrand.l2[l] = dipnorm * np.trapz(integrand.d2, dx=ddip)
        if nw == 1:
            integrand.w[w] = integrand.l[0]
            integrand.w2[w] = integrand.l2[0]
        else:
            integrand.w[w] = np.trapz(peps*integrand.l, dx=d_eps)
            integrand.w2[w] = np.trapz(peps*integrand.l2, dx=d_eps)
    if nw == 1:
        avg = integrand.w[0]
        var = integrand.w2[0] - avg**2
    else:
        avg = np.trapz(peps*integrand.w, dx=d_eps)
        var = np.trapz(peps*integrand.w2, dx=d_eps) - avg**2

    return var, avg

def main():

    start_datetime = datetime.datetime.now().isoformat()

    config_file = sys.argv[-1]
    if not os.path.exists(config_file):
        raise Exception("Config file %s doesn't exist" % (config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)
    section = Config.sections()[0]

    NP = Config.getint(section, 'NP')
    filebase = Config.get(section, 'filebase')
    try:
        datadir = Config.get(section, 'datadir')
    except:
        datadir = 'data'
    rup_dim_model = Config.get(section, 'rup_dim_model')
    mech = Config.get(section, 'mech')
    LW = Config.getboolean(section, 'LW')
    AR = Config.getfloat(section, 'AR')
    ndip = Config.getint(section, 'ndip')
    mindip_deg = Config.getfloat(section, 'mindip_deg')
    maxdip_deg = Config.getfloat(section, 'maxdip_deg')
    ntheta = Config.getint(section, 'ntheta')
    nxny = Config.getint(section, 'nxny')
    nz = Config.getint(section, 'nz')
    neps = Config.getint(section, 'neps')
    trunc = Config.getfloat(section, 'trunc')
    minmag = Config.getfloat(section, 'minmag')
    maxmag = Config.getfloat(section, 'maxmag')
    dmag = Config.getfloat(section, 'dmag')
    minepi = Config.getfloat(section, 'minepi')
    maxepi = Config.getfloat(section, 'maxepi')
    nepi = Config.getint(section, 'nepi')
    min_seis_depth = Config.getfloat(section, 'min_seis_depth')
    max_seis_depth = Config.getfloat(section, 'max_seis_depth')

    if os.path.isdir(datadir) is False:
        os.mkdir(datadir)

    if LW is True:
        filename = "%s_%s_mech%s_LW_seis%g_%g" % \
                    (filebase, rup_dim_model, mech, min_seis_depth, max_seis_depth)
    else:
        filename = "%s_%s_mech%s_ar%.1f_seis%g_%g" % \
                    (filebase, rup_dim_model, mech, AR, min_seis_depth, max_seis_depth)
    filename = filename.replace(".", "p")
    filename = os.path.join(datadir, filename)

    mindip = mindip_deg * np.pi / 180.0
    maxdip = maxdip_deg * np.pi / 180.0
    M = np.arange(minmag, maxmag + 0.001, dmag)
    Repi = np.logspace(np.log10(minepi), np.log10(maxepi), nepi)
    #Repi = np.hstack((np.logspace(np.log10(0.1), np.log10(1), 6)[0:-1],
    #                  np.logspace(np.log10(1), np.log10(1000), 31)))

    for iP in range(0, NP):
        if os.fork() == 0:
            MeanVar_MR(rup_dim_model=rup_dim_model, mech=mech, LW=LW, AR=AR,
                       ndip=ndip, mindip=mindip, maxdip=maxdip,
                       min_seis_depth=min_seis_depth, max_seis_depth=max_seis_depth,
                       ntheta=ntheta, nxny=nxny, nz=nz, neps=neps, trunc=trunc,
                       NP=NP, iP=iP, filename=filename,
                       M=M, Repi=Repi)

    for iP in range(0, NP):
        pid, status = os.waitpid(-1, 0)


    fr = open('%s_Hw.csv' % (filename), 'w')
    fv = open('%s_HwVar.csv' % (filename), 'w')

    fr.write('# Program: %s\n' % sys.argv[0])
    fv.write('# Program: %s\n' % sys.argv[0])
    fr.write('# Config file: %s\n' % config_file)
    fv.write('# Config file: %s\n' % config_file)
    fr.write('# Process start: %s\n' % start_datetime)
    fv.write('# Process start: %s\n' % start_datetime)
    fr.write('# Process finish: %s\n' % datetime.datetime.now().isoformat())
    fv.write('# Process finish: %s\n' % datetime.datetime.now().isoformat())
    fr.write('# rup_dim_model = %s, mech = %s, LW = %s, AR = %s, ndip = %d, mindip = %f, '
             'maxdip = %f, ntheta = %d, nxny = %d, nz = %d\n' \
             % (rup_dim_model, mech, LW, AR, ndip, mindip, maxdip, ntheta, nxny, nz))
    fv.write('# rup_dim_model = %s, mech = %s, LW = %s, AR = %s, ndip = %d, mindip = %f, '
             'maxdip = %f, ntheta = %d, nxny = %d, nz = %d\n' \
             % (rup_dim_model, mech, LW, AR, ndip, mindip, maxdip, ntheta, nxny, nz))
    fr.write('# neps = %d, trunc = %f, min_seis_depth = %f, max_seis_depth = %f\n'
             % (neps, trunc, min_seis_depth, max_seis_depth))
    fv.write('# neps = %d, trunc = %f, min_seis_depth = %f, max_seis_depth = %f\n' \
             % (neps, trunc, min_seis_depth, max_seis_depth))

    fr.write('"Repi_km",')
    fv.write('"Repi_km",')

    nmag = np.size(M)
    for j in range(0, nmag):
        fr.write('"R%g"' % M[j])
        fv.write('"V%g"' % M[j])
        if j < nmag - 1:
            fr.write(',')
            fv.write(',')
    fr.write("\n")
    fv.write("\n")

    rat_file = [None] * NP
    var_file = [None] * NP
    frs = [None] * NP
    fvs = [None] * NP
    for iP in range(0, NP):
        rat_file[iP] = '%sHw_%02d.csv' % (filename, iP)
        var_file[iP] = '%sHwVar_%02d.csv' % (filename, iP)
        frs[iP] = open(rat_file[iP], 'r')
        fvs[iP] = open(var_file[iP], 'r')

    for line in frs[0]:
        fr.write(line)
        line = fvs[0].readline()
        fv.write(line)
        for iP in range(1, NP):
            line = frs[iP].readline()
            if line:
                fr.write(line)
            line = fvs[iP].readline()
            if line:
                fv.write(line)

    for iP in range(0, NP):
        frs[iP].close()
        fvs[iP].close()
        os.unlink(rat_file[iP])
        os.unlink(var_file[iP])

    fr.close()
    fv.close()

if __name__ == "__main__":

    start_time = time.time()

    main()

    print("Total execution time %f seconds." % (time.time() - start_time))
    sys.exit(0)

