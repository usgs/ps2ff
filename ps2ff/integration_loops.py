#!/usr/bin/env python3


from __future__ import division
from __future__ import print_function

import sys
import datetime

import numpy as np

from ps2ff.magnitude_scaling import dimensions_from_magnitude
from ps2ff.magnitude_scaling import compute_epsilon


def mag_dist_loop(what, ndip=19, mindip=0.0, maxdip=np.pi/2.0,
                  min_seis_depth=0, max_seis_depth=20,
                  rup_dim_model='WC94', mech="A", LW=True,
                  AR=1, ntheta=73, nxny=100, neps=10, trunc=3,
                  NP=1, iP=0, filename='Repi-to-Rjb-',
                  M=6.0, Repi=100, nz=20):
    """
    This function loops over M and R, calling rjb_inner_loop
    and writes out progress information.

    Args:
        what (str): What to ingegrate? Currently supported values are
            - 'Rrup'
            - 'Rjb'
        ndip (int): Number of integration steps for dip.
        mindip (float): The minimum rupture dip in degrees (0 to 90).
        maxdip (float): The maximum rupture dip in degrees (0 to 90).
        min_seis_depth (float): The minimum seismogenic depth (km).
        max_seis_depth (float): The maximum seismogenic depth (km).
        rup_dim_model (str): String indicating the model for compputing the
            rupture dimensions from magnitude. Supported values are:
                - 'WC94'
                - 'S14'
                - 'HB08'
        mech (str): Optional string indicating earthquake mechanism, used by
            some of the models. Anything other than 'R', 'N', 'SS', or 'A'
            (the default).
        LW (bool): Compute length and width from magnitude, and integrate
            across them individually. Alternative is to assume an aspect
            ratio.
        AR (float): Aspect ratio (length/width).
        ntheta (int): Number of integration steps for theta.
        nxny (int): Number of integration steps in the x and y direction.
        neps (int): Number of integration steps for epsilon.
        trunc (float): Epsilon truncation level.
        NP (int): Number of forked processes.
        iP (int): Multiple process index.
        filename (str): Output file name.
        M (float): Earthquake magnitude.
        Repi (float): Epicentral distance (km).
        nz (int): Number of integration steps in depth. Only used
            for Rrup calculations.
    """
    if what not in ['Rjb', 'Rrup']:
        raise Exception('Unsupported distance type.')

    nepi = np.size(Repi)

    if NP != 1:
        ii = range(iP, nepi, NP)
        Repi = Repi[ii]
        nepi = np.size(Repi)

    nmag = np.size(M)
    ratio = np.zeros((nepi, nmag))
    variance = np.zeros((nepi, nmag))
    for i in range(nepi):
        for j in range(nmag):
            if what == 'Rjb':
                dist_var, dist_avg = rjb_inner_loop(
                    M[j], Repi[i], ndip=ndip, mindip=mindip, maxdip=maxdip,
                    min_seis_depth=min_seis_depth,
                    max_seis_depth=max_seis_depth,
                    rup_dim_model=rup_dim_model, mech=mech,
                    LW=LW, AR=AR, ntheta=ntheta, nxny=nxny, neps=neps,
                    trunc=trunc)
            elif what == 'Rrup':
                dist_var, dist_avg = rrup_inner_loop(
                    M[j], Repi[i], ndip=ndip, mindip=mindip, maxdip=maxdip,
                    min_seis_depth=min_seis_depth,
                    max_seis_depth=max_seis_depth,
                    rup_dim_model=rup_dim_model, mech=mech,
                    LW=LW, AR=AR, ntheta=ntheta, nxny=nxny, nz=nz,
                    neps=neps,
                    trunc=trunc)

            ratio[i, j] = dist_avg / Repi[i]
            variance[i, j] = dist_var
            print("        Proc %d j=%d of %d %s" % (iP, j+1, nmag,
                  datetime.datetime.now().isoformat()))
        print("Proc %d done with %d of %d distances %s" % (iP, i+1, nepi,
              datetime.datetime.now().isoformat()))

    fr = open('%sRatios_%02d.csv' % (filename, iP), 'w')
    fv = open('%sVar_%02d.csv' % (filename, iP), 'w')

    for i in range(nepi):
        fr.write("%f," % Repi[i])
        fv.write("%f," % Repi[i])
        for j in range(nmag):
            fr.write("%f" % ratio[i, j])
            fv.write("%f" % variance[i, j])
            if j < nmag - 1:
                fr.write(",")
                fv.write(",")
        fr.write("\n")
        fv.write("\n")

    fr.close()
    fv.close()
    sys.exit(0)


def rjb_inner_loop(M, Repi, ndip=19, mindip=0.0, maxdip=np.pi/2.0,
                   min_seis_depth=0,  max_seis_depth=20,
                   rup_dim_model='WC94', mech="A", LW=True, AR=1,
                   ntheta=73, nxny=50,
                   neps=10, trunc=3):
    """
    This function evaluates the Rjb mean and var integral
    for a single M/R pair, looping over:
       - dip
       - dx, dy (location of hypocenter on fault)
       - theta (angle to fault)
       - epsilon (dummy variable for L/W/A integration)
    We do this so that parallizaiton is simple: this function can be forked
    onto different cores.

    See the README for argument definitions.
    """
    max_seis_thickness = max_seis_depth - min_seis_depth
    # Same as before, but position of P is random so have to integrate from
    # theta = 0 to 2pi and then for dx = 0 to L and dy = 0 to SW
    #  and dip from 0 to pi/2
    if ndip != 1:
        dip = np.linspace(mindip, maxdip, ndip)
        ddip = dip[1] - dip[0]
        dipnorm = 1.0 / (maxdip - mindip)
    else:
        dip = np.array([mindip])

    length, sig_length, width, sigw, area, sig_area = \
        dimensions_from_magnitude(M, rup_dim_model, neps, trunc, mech)

    if LW is True:
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

        width_matrix = np.tile(np.sqrt(area/AR), (ndip, 1))
        sindip = np.tile(np.sin(dip).reshape(-1, 1), (1, neps))
        rup_z = width_matrix * sindip
        indxx = rup_z > max_seis_thickness
        width_matrix[indxx] = max_seis_thickness / sindip[indxx]
        length_matrix = np.tile(area, (ndip, 1)) / width_matrix
        nl = 1
        nw = neps

    theta = np.linspace(0, 2*np.pi, ntheta)  # in rad
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

    epsmid, peps, d_eps = compute_epsilon(neps, trunc)

    integrand_width = np.zeros(nw) + np.nan
    integrand_width2 = np.zeros(nw) + np.nan
    integrand_length = np.zeros(nl)
    integrand_length2 = np.zeros(nl)
    integrand_dip = np.zeros(ndip)
    integrand_dip2 = np.zeros(ndip)
    Rjb = np.zeros(ntheta)

    for w in range(nw):
        for l in range(nl):
            for k in range(ndip):
                if LW is True:
                    ll = l
                else:
                    width = width_matrix[k, :]
                    length = length_matrix[k, :]
                    ll = w
                one_over_ll = 1.0 / length[ll]

                ny = nxny
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
                    nx = nxny
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
            if ndip == 1:
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


def rrup_inner_loop(M, Repi, ndip=19, mindip=0.0, maxdip=np.pi/2.0,
                    min_seis_depth=0,  max_seis_depth=20,
                    rup_dim_model='WC94', mech="A", LW=True, AR=1,
                    ntheta=73, nxny=50, nz=20,
                    neps=10, trunc=3):
    """
    """
    max_seis_thickness = max_seis_depth - min_seis_depth
    # Same as before, but position of P is random so have to integrate from
    # theta = 0 to 2pi and then for dx = 0 to L and dy = 0 to SW
    # and dip from 0 to pi/2
    if ndip != 1:
        dip = np.linspace(mindip, maxdip, ndip)
        ddip = dip[1] - dip[0]
        dipnorm = 1.0 / (maxdip - mindip)
    else:
        dip = np.array([mindip])

    length, sig_length, width, sigw, area, sig_area = \
        dimensions_from_magnitude(M, rup_dim_model, neps, trunc, mech)

    if LW is True:
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

        width_matrix = np.tile(np.sqrt(area/AR), (ndip, 1))
        sindip = np.tile(np.sin(dip).reshape(-1, 1), (1, neps))
        rup_z = width_matrix * sindip
        indxx = rup_z > max_seis_thickness
        width_matrix[indxx] = max_seis_thickness / sindip[indxx]
        length_matrix = np.tile(area, (ndip, 1)) / width_matrix
        nl = 1
        nw = neps

    theta = np.linspace(0, 2*np.pi, ntheta)  # in rad
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

    epsmid, peps, d_eps = compute_epsilon(neps, trunc)

    integrand_width = np.zeros(nw) + np.nan
    integrand_width2 = np.zeros(nw) + np.nan
    integrand_length = np.zeros(nl)
    integrand_length2 = np.zeros(nl)
    integrand_dip = np.zeros(ndip)
    integrand_dip2 = np.zeros(ndip)
    integrand_depth = np.zeros(nz)
    integrand_depth2 = np.zeros(nz)
    Rjb = np.zeros(ntheta)
    Rrup = np.zeros(ntheta)
    Rrupp = np.zeros(ntheta)
    Ry = np.zeros(ntheta)

    for w in range(nw):
        for l in range(nl):
            for k in range(ndip):
                if LW is True:
                    ll = l
                else:
                    width = width_matrix[k, :]
                    length = length_matrix[k, :]
                    ll = w
                one_over_ll = 1.0 / length[ll]

                # Since we still use linspace, dy won't be exactly minx.
                # Also put a hard bound on minimum ny to be 2 so that dy
                # makes sense.
                ny = nxny
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
                    nx = nxny
                    x = np.linspace(0, SW, nx)
                    dx = x[1] - x[0]
                # Calclate range of Ztor
                ZtorMax = max_seis_depth - width[w]*np.sin(dip[k])
                if np.allclose(ZtorMax, 0) or nz == 1:
                    nz = 1
                    dz = 0
                    Ztor = np.linspace(0, ZtorMax, nz)
                else:
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
            if ndip == 1:
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
