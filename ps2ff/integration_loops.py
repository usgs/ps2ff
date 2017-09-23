#!/usr/bin/env python3


from __future__ import division
from __future__ import print_function

import sys
import datetime

import numpy as np

from ps2ff.magnitude_scaling import dimensions_from_magnitude
from ps2ff.magnitude_scaling import compute_epsilon
from ps2ff.constants import DistType, Mechanism, MagScaling


def mag_dist_loop(conf, iP=None, filename=None, M=None, Repi=None):
    """
    This function loops over M and R, calling rjb_inner_loop
    and writes out progress information.

    Args:
        conf(ConfigObj): The configuration info. See
            `ps2ff/data/configspec.conf`.
        iP (int): Multiple process index.
        filename (str): Output file name.
        M (numpy.ndarray): Earthquake magnitudes.
        Repi (numpy.ndarray): Epicentral distances (km).
    """
    nepi = np.size(Repi)

    if conf['NP'] != 1:
        ii = range(iP, nepi, conf['NP'])
        Repi = Repi[ii]
        nepi = np.size(Repi)

    nmag = np.size(M)
    ratio = np.zeros((nepi, nmag))
    variance = np.zeros((nepi, nmag))
    for i in range(nepi):
        for j in range(nmag):
            if conf['what'] is DistType.Rjb:
                dist_var, dist_avg = rjb_inner_loop(M[j], Repi[i], conf)
            elif conf['what'] is DistType.Rrup:
                dist_var, dist_avg = rrup_inner_loop(M[j], Repi[i], conf)
            else:
                raise TypeError('Unknown distance type: %s' % (conf['what']))

            ratio[i, j] = dist_avg / Repi[i]
            variance[i, j] = dist_var
            print("        Proc %d j=%d of %d %s" % (iP, j+1, nmag,
                  datetime.datetime.now().isoformat()))
        print("Proc %d done with %d of %d distances %s" % (iP, i+1, nepi,
              datetime.datetime.now().isoformat()))

    for fname, farr in (('Ratios', ratio), ('Var', variance)):
        f = open('%s%s_%02d.csv' % (filename, fname, iP), 'w')
        for i in range(nepi):
            f.write("%f," % Repi[i])
            for j in range(nmag):
                f.write("%f" % farr[i, j])
                if j < nmag - 1:
                    f.write(",")
            f.write("\n")
        f.close()
    sys.exit(0)


def single_event_loop(conf, iP=0, rjb_filename='junk',
                      rrup_filename='junk', M=6.0, Repi=None):
    """
    Args:
        conf(ConfigObj): The configuration info. See
            `ps2ff/data/configspec.conf`.
        iP (int): Multiple process index.
        rjb_filename (str): Output file name for Rjb results.
        rrup_filename (str): Output file name for Rrup results.
        M (float): Earthquake magnitude.
        Repi (float): Epicentral distance (km).
    """
    nepi = np.size(Repi)

    if conf['NP'] != 1:
        ii = range(iP, conf['nepi'], conf['NP'])
        Repi = Repi[ii]
        nepi = np.size(Repi)
    print('ip=%d nepi=%d Repi: ' % (iP, nepi), Repi)

    if conf['bytheta'] is True:
        Rrup_ratio = np.zeros((nepi, conf['ntheta']))
        Rrup_variance = np.zeros((nepi, conf['ntheta']))
        Rjb_ratio = np.zeros((nepi, conf['ntheta']))
        Rjb_variance = np.zeros((nepi, conf['ntheta']))
    else:
        Rrup_ratio = np.zeros((nepi, 1))
        Rrup_variance = np.zeros((nepi, 1))
        Rjb_ratio = np.zeros((nepi, 1))
        Rjb_variance = np.zeros((nepi, 1))

    for i in range(0, nepi):
        if conf['bytheta'] is True:
            theta = np.linspace(0, np.pi*2, conf['ntheta'])  # in rad
            for j in range(0, conf['ntheta']):
                rrup_var, rrup_avg, rjb_var, rjb_avg = \
                    single_event_inner_loop(conf, Repi[i],
                                            theta=np.array([theta[j]]),
                                            ntheta=1)
                Rrup_ratio[i, j] = rrup_avg / Repi[i]
                Rrup_variance[i, j] = rrup_var
                Rjb_ratio[i, j] = rjb_avg / Repi[i]
                Rjb_variance[i, j] = rjb_var
                print("        Proc %d j=%d of %d %s" % (iP, j+1,
                      conf['ntheta'], datetime.datetime.now().isoformat()))
                print("Proc %d done with %d of %d distances %s" %
                      (iP, i+1, nepi, datetime.datetime.now().isoformat()))
        else:
            rrup_var, rrup_avg, rjb_var, rjb_avg = \
                single_event_inner_loop(conf, Repi[i], theta=0,
                                        ntheta=conf['ntheta'])
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
        if conf['bytheta'] is True:
            for j in range(0, conf['ntheta']):
                fr_rrup.write("%f" % Rrup_ratio[i, j])
                fv_rrup.write("%f" % Rrup_variance[i, j])
                fr_rjb.write("%f" % Rjb_ratio[i, j])
                fv_rjb.write("%f" % Rjb_variance[i, j])
                if j < conf['ntheta'] - 1:
                    fr_rrup.write(",")
                    fv_rrup.write(",")
                    fr_rjb.write(",")
                    fv_rjb.write(",")
            fr_rrup.write("\n")
            fv_rrup.write("\n")
            fr_rjb.write("\n")
            fv_rjb.write("\n")
        else:
            fr_rrup.write("%f\n" % Rrup_ratio[i, 0])
            fv_rrup.write("%f\n" % Rrup_variance[i, 0])
            fr_rjb.write("%f\n" % Rjb_ratio[i, 0])
            fv_rjb.write("%f\n" % Rjb_variance[i, 0])

    fr_rrup.close()
    fv_rrup.close()
    fr_rjb.close()
    fv_rjb.close()
    sys.exit(0)


def rjb_inner_loop(M, Repi, conf):
    """
    This function evaluates the Rjb mean and var integral
    for a single M/R pair, looping over:
       - dip
       - dx, dy (location of hypocenter on fault)
       - theta (angle to fault)
       - epsilon (dummy variable for L/W/A integration)
    We do this so that parallizaiton is simple: this function can be forked
    onto different cores.

    Args:
        M (float): Earthquake magnitude.
        Repi (float): Epicentral distance (km).
        conf(ConfigObj): The configuration info. See
            `ps2ff/data/configspec.conf`.

    Returns:
        tuple: Rjb variance, mean Rjb.
    """
    max_seis_thickness = conf['max_seis_depth'] - conf['min_seis_depth']
    if conf['ndip'] != 1:
        dip = np.linspace(conf['mindip'], conf['maxdip'], conf['ndip'])
        ddip = dip[1] - dip[0]
        dipnorm = 1.0 / (conf['maxdip'] - conf['mindip'])
    else:
        dip = np.array([conf['mindip']])

    length, sig_length, width, sigw, area, sig_area = \
        dimensions_from_magnitude(M, conf['rup_dim_model'], conf['neps'],
                                  conf['trunc'], conf['mech'])

    if conf['LW'] is True:
        if length is None or width is None:
            raise Exception('Selected model does not support direct '
                            'estimation of length and width. Either '
                            'select an one that does, or use an assumed '
                            'aspect ratio.')
        else:
            nl = len(length)
            nw = len(width)
    else:
        # Trick the L and W loops to handle a one to one mapping
        # between L and W, and each L/W pair corresponds to a
        # single M looping over epsilon; specify length and
        # width constrained by seismogenic depth.

        width_matrix = np.tile(np.sqrt(area/conf['AR']), (conf['ndip'], 1))
        sindip = np.tile(np.sin(dip).reshape(-1, 1), (1, conf['neps']))
        rup_z = width_matrix * sindip
        indxx = rup_z > max_seis_thickness
        width_matrix[indxx] = max_seis_thickness / sindip[indxx]
        length_matrix = np.tile(area, (conf['ndip'], 1)) / width_matrix
        nl = 1
        nw = conf['neps']

    theta = np.linspace(0, 2*np.pi, conf['ntheta'])  # in rad
    dt = theta[1] - theta[0]
    # origin defined at o:
    #
    #  + +------+
    #  | |      |
    #  L |      |
    #  | |      |
    #  + o------+
    #    +--SW--+
    #
    one_over_2pi = 1.0/2.0/np.pi

    epsmid, peps, d_eps = compute_epsilon(conf['neps'], conf['trunc'])

    integrand_width = np.zeros(nw) + np.nan
    integrand_width2 = np.zeros(nw) + np.nan
    integrand_length = np.zeros(nl)
    integrand_length2 = np.zeros(nl)
    integrand_dip = np.zeros(conf['ndip'])
    integrand_dip2 = np.zeros(conf['ndip'])
    Rjb = np.zeros(conf['ntheta'])

    for w in range(nw):
        for l in range(nl):
            for k in range(conf['ndip']):
                if conf['LW'] is True:
                    ll = l
                else:
                    width = width_matrix[k, :]
                    length = length_matrix[k, :]
                    ll = w
                one_over_ll = 1.0 / length[ll]

                ny = conf['nxny']
                y = np.linspace(0, length[ll], ny)
                dy = y[1] - y[0]

                integrand_y = np.zeros(ny)
                integrand_y2 = np.zeros(ny)

                # Fault width projected to surface:
                if np.allclose(np.cos(dip[k]), 0):
                    SW = 0
                    nx = 1
                    x = np.linspace(0, SW, nx)
                    dx = 0
                else:
                    SW = width[w] * np.cos(dip[k])
                    one_over_sw = 1.0 / SW

                    # Since we still use linspace, dx won't be exactly minx.
                    # Also put a hard bound on minimum nx to be 2 so that dx
                    # makes sense.
                    nx = conf['nxny']
                    x = np.linspace(0, SW, nx)
                    dx = x[1] - x[0]

                integrand_x = np.zeros(nx)
                integrand_x2 = np.zeros(nx)
                for i in range(nx):
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
                    for j in range(ny):
                        yi = y[j] + Repi * np.sin(theta)

                        cca = yi > length[ll]
                        ccb = (yi >= 0) & (yi <= length[ll])
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
                        # The original version below. The above two chunks
                        # are an optimization of the following:
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
                        yy = yi[c1] - length[ll]
                        Rjb[c1] = np.sqrt(xx*xx + yy*yy)
                        Rjb[c2] = yi[c2] - length[ll]
                        xx = xj[c3] - SW
                        yy = yi[c3] - length[ll]
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
                        Rjb2 = Rjb * Rjb
                        integrand_y[j] = one_over_2pi * np.trapz(Rjb, dx=dt)
                        integrand_y2[j] = one_over_2pi * np.trapz(Rjb2, dx=dt)
                    integrand_x[i] = one_over_ll * np.trapz(integrand_y, dx=dy)
                    integrand_x2[i] = one_over_ll * \
                        np.trapz(integrand_y2, dx=dy)
                if dx == 0:
                    integrand_dip[k] = integrand_x[0]
                    integrand_dip2[k] = integrand_x2[0]
                else:
                    integrand_dip[k] = one_over_sw * \
                        np.trapz(integrand_x, dx=dx)
                    integrand_dip2[k] = one_over_sw * \
                        np.trapz(integrand_x2, dx=dx)
            if conf['ndip'] == 1:
                integrand_length[l] = integrand_dip[0]
                integrand_length2[l] = integrand_dip2[0]
            else:
                integrand_length[l] = dipnorm * \
                    np.trapz(integrand_dip, dx=ddip)
                integrand_length2[l] = dipnorm * \
                    np.trapz(integrand_dip2, dx=ddip)
        if nw == 1:
            integrand_width[w] = integrand_length[0]
            integrand_width2[w] = integrand_length2[0]
        else:
            integrand_width[w] = np.trapz(peps*integrand_length, dx=d_eps)
            integrand_width2[w] = np.trapz(peps*integrand_length2, dx=d_eps)
    if nw == 1:
        Rjb_avg = integrand_width[0]
        Rjb_var = integrand_width2[0] - Rjb_avg**2
    else:
        Rjb_avg = np.trapz(peps*integrand_width, dx=d_eps)
        Rjb_var = np.trapz(peps*integrand_width2, dx=d_eps) - Rjb_avg**2

    return Rjb_var, Rjb_avg


def rrup_inner_loop(M, Repi, conf):
    """
    This function evaluates the Rrup mean and var integral
    for a single M/R pair, looping over:
       - dip
       - dx, dy (location of hypocenter on fault)
       - theta (angle to fault)
       - epsilon (dummy variable for L/W/A integration)
    We do this so that parallizaiton is simple: this function can be forked
    onto different cores.

    Args:
        M (float): Earthquake magnitude.
        Repi (float): Epicentral distance (km).
        conf(ConfigObj): The configuration info. See
            `ps2ff/data/configspec.conf`.

    Returns:
        tuple: Rjb variance, mean Rjb.
    """
    max_seis_thickness = conf['max_seis_depth'] - conf['min_seis_depth']
    if conf['ndip'] != 1:
        dip = np.linspace(conf['mindip'], conf['maxdip'], conf['ndip'])
        ddip = dip[1] - dip[0]
        dipnorm = 1.0 / (conf['maxdip'] - conf['mindip'])
    else:
        dip = np.array([conf['mindip']])

    length, sig_length, width, sigw, area, sig_area = \
        dimensions_from_magnitude(M, conf['rup_dim_model'], conf['neps'],
                                  conf['trunc'], conf['mech'])

    if conf['LW'] is True:
        if length is None or width is None:
            raise Exception('Selected model does not support direct '
                            'estimation of length and width. Either '
                            'select an one that does, or use an assumed '
                            'aspect ratio.')
        else:
            nl = len(length)
            nw = len(width)
    else:
        # Trick the L and W loops to handle a one to one mapping
        # between L and W, and each L/W pair corresponds to a
        # single M looping over epsilon; specify length and
        # width constrained by seismogenic depth.

        width_matrix = np.tile(np.sqrt(area/conf['AR']), (conf['ndip'], 1))
        sindip = np.tile(np.sin(dip).reshape(-1, 1), (1, conf['neps']))
        rup_z = width_matrix * sindip
        indxx = rup_z > max_seis_thickness
        width_matrix[indxx] = max_seis_thickness / sindip[indxx]
        length_matrix = np.tile(area, (conf['ndip'], 1)) / width_matrix
        nl = 1
        nw = conf['neps']

    theta = np.linspace(0, 2*np.pi, conf['ntheta'])  # in rad
    dt = theta[1] - theta[0]

    # origin defined at o:
    #
    #  + +------+
    #  | |      |
    #  L |      |
    #  | |      |
    #  + o------+
    #    +--SW--+
    #
    one_over_2pi = 1.0/2.0/np.pi

    epsmid, peps, d_eps = compute_epsilon(conf['neps'], conf['trunc'])

    integrand_width = np.zeros(nw) + np.nan
    integrand_width2 = np.zeros(nw) + np.nan
    integrand_length = np.zeros(nl)
    integrand_length2 = np.zeros(nl)
    integrand_dip = np.zeros(conf['ndip'])
    integrand_dip2 = np.zeros(conf['ndip'])
    integrand_depth = np.zeros(conf['nz'])
    integrand_depth2 = np.zeros(conf['nz'])
    Rjb = np.zeros(conf['ntheta'])
    Rrup = np.zeros(conf['ntheta'])
    Rrupp = np.zeros(conf['ntheta'])
    Ry = np.zeros(conf['ntheta'])

    for w in range(nw):
        for l in range(nl):
            for k in range(conf['ndip']):
                if conf['LW'] is True:
                    ll = l
                else:
                    width = width_matrix[k, :]
                    length = length_matrix[k, :]
                    ll = w
                one_over_ll = 1.0 / length[ll]

                # Since we still use linspace, dy won't be exactly minx.
                # Also put a hard bound on minimum ny to be 2 so that dy
                # makes sense.
                ny = conf['nxny']
                y = np.linspace(0, length[ll], ny)
                dy = y[1] - y[0]

                integrand_y = np.zeros(ny)
                integrand_y2 = np.zeros(ny)

                # Fault width projected to surface:
                if np.allclose(np.cos(dip[k]), 0):
                    SW = 0
                    nx = 1
                    x = np.linspace(0, SW, nx)
                    dx = 0
                else:
                    SW = width[w] * np.cos(dip[k])
                    one_over_sw = 1.0 / SW

                    # Since we still use linspace, dx won't be exactly minx.
                    # Also put a hard bound on minimum nx to be 2 so that dx
                    # makes sense.
                    nx = conf['nxny']
                    x = np.linspace(0, SW, nx)
                    dx = x[1] - x[0]
                # Calclate range of Ztor
                ZtorMax = conf['max_seis_depth'] - width[w]*np.sin(dip[k])
                if np.allclose(ZtorMax, 0) or conf['nz'] == 1: # Should 0 be min_seis_depth
                    nz = 1
                    dz = 0
                    Ztor = np.array([0])        # should be min_seis_depth ???
                else:
                    nz = conf['nz']
                    one_over_ztormax = 1.0/ZtorMax
                    Ztor = np.linspace(0, ZtorMax, nz)
                    dz = Ztor[1]-Ztor[0]

                integrand_x = np.zeros(nx)
                integrand_x2 = np.zeros(nx)
                for z in range(nz):
                    for i in range(nx):
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
                        for j in range(ny):
                            yi = y[j] + Repi * np.sin(theta)

                            cca = yi > length[ll]
                            ccb = (yi >= 0) & (yi <= length[ll])
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
                            yy = yi[c1] - length[ll]
                            Rjb[c1] = np.sqrt(xx*xx + yy*yy)
                            Rjb[c2] = yi[c2] - length[ll]
                            xx = xj[c3] - SW
                            yy = yi[c3] - length[ll]
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
                            Ry[c1 | c2 | c3] = yi[c1 | c2 | c3] - length[ll]
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
                                r2 = (tmp <= Rx) & \
                                    (Rx <= (tmp + width[w]/np.cos(dip[k])))
                                Rrupp[r2] = Rx[r2]*np.sin(dip[k]) + \
                                    Ztor[z]*np.cos(dip[k])
                                r3 = Rx > (tmp + width[w]/np.cos(dip[k]))
                                Rrupp[r3] = np.sqrt((Rx[r3] -
                                    width[w]*np.cos(dip[k]))**2 + (Ztor[z] +
                                    width[w]*np.sin(dip[k]))**2)

                            Rrup = np.sqrt(Rrupp**2 + Ry**2)
                            Rrup2 = Rrup * Rrup

                            integrand_y[j] = one_over_2pi * \
                                np.trapz(Rrup, dx=dt)
                            integrand_y2[j] = one_over_2pi * \
                                np.trapz(Rrup2, dx=dt)

                        integrand_x[i] = one_over_ll * \
                            np.trapz(integrand_y, dx=dy)
                        integrand_x2[i] = one_over_ll * \
                            np.trapz(integrand_y2, dx=dy)

                    if dx == 0:
                        integrand_depth[z] = integrand_x[0]
                        integrand_depth2[z] = integrand_x2[0]
                    else:
                        integrand_depth[z] = one_over_sw * \
                            np.trapz(integrand_x, dx=dx)
                        integrand_depth2[z] = one_over_sw * \
                            np.trapz(integrand_x2, dx=dx)
                # end z loop
                if dz == 0:
                    integrand_dip[k] = integrand_depth[0]
                    integrand_dip2[k] = integrand_depth2[0]
                else:
                    integrand_dip[k] = one_over_ztormax * \
                        np.trapz(integrand_depth, dx=dz)
                    integrand_dip2[k] = one_over_ztormax * \
                        np.trapz(integrand_depth2, dx=dz)
            if conf['ndip'] == 1:
                integrand_length[l] = integrand_dip[0]
                integrand_length2[l] = integrand_dip2[0]
            else:
                integrand_length[l] = dipnorm * \
                    np.trapz(integrand_dip, dx=ddip)
                integrand_length2[l] = dipnorm * \
                    np.trapz(integrand_dip2, dx=ddip)
        if nw == 1:
            integrand_width[w] = integrand_length[0]
            integrand_width2[w] = integrand_length2[0]
        else:
            integrand_width[w] = np.trapz(peps*integrand_length, dx=d_eps)
            integrand_width2[w] = np.trapz(peps*integrand_length2, dx=d_eps)
    if nw == 1:
        Rrup_avg = integrand_width[0]
        Rrup_var = integrand_width2[0] - Rrup_avg**2
    else:
        Rrup_avg = np.trapz(peps*integrand_width, dx=d_eps)
        Rrup_var = np.trapz(peps*integrand_width2, dx=d_eps) - Rrup_avg**2

    return Rrup_var, Rrup_avg


def single_event_inner_loop(conf, Repi, theta=0, ntheta=73):
    """
    Args:
        conf(ConfigObj): The configuration info. See
           `ps2ff/data/configspec.conf`.
        Repi (float): Epicentral distance (km).
        theta (float): Source-to-site angle (radians).
        ntheta (int): Number of integration steps for theta; used if `bytheta`
            is True.

    Returns:
        tuple: Rrup variance, Rrup mean, Rjb variance, Rjb mean
    """
    if conf['ndip'] != 1:
        dip = np.linspace(conf['mindip'], conf['maxdip'], conf['ndip'])
        ddip = dip[1] - dip[0]
        dipnorm = 1.0 / (conf['maxdip'] - conf['mindip'])
    else:
        dip = np.array([conf['mindip']])

    length, sig_length, width, sigw, area, sig_area = \
        dimensions_from_magnitude(conf['M'], conf['rup_dim_model'],
                                  conf['neps'], conf['trunc'], conf['mech'])
#    area = area[0]  # fix dimensions

    if conf['bytheta'] is False:
        theta = np.linspace(0, 2*np.pi, ntheta)  # in rad
        dt = theta[1] - theta[0]

    one_over_2pi = 1.0/2.0/np.pi

    epsmid, peps, d_eps = compute_epsilon(conf['neps'], conf['trunc'])

    RrupIntegrand_a = np.zeros(conf['neps']) + np.nan
    RrupIntegrand_a2 = np.zeros(conf['neps']) + np.nan
    RrupIntegrand_d = np.zeros(conf['ndip'])
    RrupIntegrand_d2 = np.zeros(conf['ndip'])
    RrupIntegrand_y = np.zeros(conf['nxny'])
    RrupIntegrand_y2 = np.zeros(conf['nxny'])
    RrupIntegrand_x = np.zeros(conf['nxny'])
    RrupIntegrand_x2 = np.zeros(conf['nxny'])

    RjbIntegrand_a = np.zeros(conf['neps']) + np.nan
    RjbIntegrand_a2 = np.zeros(conf['neps']) + np.nan
    RjbIntegrand_d = np.zeros(conf['ndip'])
    RjbIntegrand_d2 = np.zeros(conf['ndip'])
    RjbIntegrand_y = np.zeros(conf['nxny'])
    RjbIntegrand_y2 = np.zeros(conf['nxny'])
    RjbIntegrand_x = np.zeros(conf['nxny'])
    RjbIntegrand_x2 = np.zeros(conf['nxny'])

    Rjb = np.zeros(ntheta)
    Rrup = np.zeros(ntheta)
    Rrupp = np.zeros(ntheta)
    Ry = np.zeros(ntheta)
    nx = conf['nxny']
    ny = conf['nxny']

    for m in range(0, conf['neps']):  # area
        W = np.sqrt(area[m]/conf['AR'])
        for k in range(0, conf['ndip']):  # dip
            if np.allclose(dip[k], 0) is False:
                ZW = W * np.sin(dip[k])  # vertical projection of W
                # Overwrite W if it extends too far and recompute SW, x, dx
                Ztor = np.max(np.array([conf['zhyp'] - ZW,
                                        conf['min_seis_depth']]))
                Zbor = np.min(np.array([conf['zhyp'] + ZW,
                                        conf['max_seis_depth']]))
                W = (Zbor - Ztor)/np.sin(dip[k])
            else:
                Ztor = conf['zhyp']

            L = area[m]/W
            SW = W * np.cos(dip[k])
            x = np.linspace(0, SW, nx)
            dx = x[1] - x[0]
            one_over_L = 1.0 / L
            one_over_sw = 1.0 / SW
            y = np.linspace(0, L, ny)
            dy = y[1] - y[0]
            for i in range(0, nx):  # x
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
                    # Compute Rx and Ry
                    Rx = xj
                    Ry[c1 | c2 | c3] = yi[c1 | c2 | c3] - L
                    Ry[c4 | c5 | c6] = 0
                    Ry[c7 | c8 | c9] = np.abs(yi[c7 | c8 | c9])
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
                        Rrupp[r3] = np.sqrt((Rx[r3] - W*np.cos(dip[k]))**2 +
                                            (Ztor + W*np.sin(dip[k]))**2)

                    Rrup = np.sqrt(Rrupp**2 + Ry**2)
                    Rrup2 = Rrup * Rrup
                    Rjb2 = Rjb * Rjb
                    if conf['bytheta'] is False:
                        RrupIntegrand_y[j] = one_over_2pi * \
                            np.trapz(Rrup, dx=dt)
                        RrupIntegrand_y2[j] = one_over_2pi * \
                            np.trapz(Rrup2, dx=dt)
                        RjbIntegrand_y[j] = one_over_2pi * \
                            np.trapz(Rjb, dx=dt)
                        RjbIntegrand_y2[j] = one_over_2pi * \
                            np.trapz(Rjb2, dx=dt)
                    else:
                        RrupIntegrand_y[j] = Rrup
                        RrupIntegrand_y2[j] = Rrup2
                        RjbIntegrand_y[j] = Rjb
                        RjbIntegrand_y2[j] = Rjb2

                RrupIntegrand_x[i] = one_over_L * \
                    np.trapz(RrupIntegrand_y, dx=dy)
                RrupIntegrand_x2[i] = one_over_L * \
                    np.trapz(RrupIntegrand_y2, dx=dy)
                RjbIntegrand_x[i] = one_over_L * \
                    np.trapz(RjbIntegrand_y, dx=dy)
                RjbIntegrand_x2[i] = one_over_L * \
                    np.trapz(RjbIntegrand_y2, dx=dy)

            RrupIntegrand_d[k] = one_over_sw * \
                np.trapz(RrupIntegrand_x, dx=dx)
            RrupIntegrand_d2[k] = one_over_sw * \
                np.trapz(RrupIntegrand_x2, dx=dx)
            RjbIntegrand_d[k] = one_over_sw * \
                np.trapz(RjbIntegrand_x, dx=dx)
            RjbIntegrand_d2[k] = one_over_sw * \
                np.trapz(RjbIntegrand_x2, dx=dx)
        if conf['ndip'] == 1:
            RrupIntegrand_a[m] = RrupIntegrand_d[0]
            RrupIntegrand_a2[m] = RrupIntegrand_d2[0]
            RjbIntegrand_a[m] = RjbIntegrand_d[0]
            RjbIntegrand_a2[m] = RjbIntegrand_d2[0]
        else:
            RrupIntegrand_a[m] = dipnorm * \
                np.trapz(RrupIntegrand_d, dx=ddip)
            RrupIntegrand_a2[m] = dipnorm * \
                np.trapz(RrupIntegrand_d2, dx=ddip)
            RjbIntegrand_a[m] = dipnorm * \
                np.trapz(RjbIntegrand_d, dx=ddip)
            RjbIntegrand_a2[m] = dipnorm * \
                np.trapz(RjbIntegrand_d2, dx=ddip)
    if conf['neps'] == 1:
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
