#!/usr/bin/env python3

import numpy as np
from scipy.stats import norm

from ps2ff.constants import MagScaling, Mechanism


def dimensions_from_magnitude(M, rup_dim_model, neps, trunc,
                              mech=Mechanism.A):
    """
    Compute dimensions of rupture from magnitude for a specified
    magnitude scaling relation.

    Args:
        M (float): Magnitude.
        rup_dim_model (MagScaling enum): Specifies the model for compputing the
            rupture dimensions from magnitude.
        neps (int): The number of steps to integrate from -trunc to +trunc.
            Larger numbers increase the accuracy of the result, but take
            longer to run.
        trunc (float): For the integration in area (or length and width), trunc
            is the truncation of the normal distribution (in units of sigma).
        mech (Mechanism enum): Optional string indicating earthquake
            mechanism, used by some of the models.

    Returns:
        tuple: A tuple containing the following, noting that some of these will
        be empty if the selected model does not provide them:

                - length: rupture length (km).
                - sig_length: standard deviation of rupture length.
                - W: rupture width (km).
                - sigw: standard devation of rupture width.
                - A: rupture area (km).
                - siga: standard deivaiton of rupture area.

    """
    epsmid, peps, d_eps = compute_epsilon(neps, trunc)
    if not isinstance(rup_dim_model, MagScaling):
        raise TypeError('rup_dim_model must be of type MagScaling')
    if not isinstance(mech, Mechanism):
        raise TypeError('mech must be of type Mechanism')

    if rup_dim_model is MagScaling.WC94:
        # Use mech to get either M-A or (M-W) and (M-R) from Wells and
        # Coppersmith.
        if mech is Mechanism.SS:
            sig_length = 0.15
            length = 10**(-2.57 + 0.62 * M + sig_length * epsmid)
            sig_width = 0.14
            width = 10**(-0.76 + 0.27 * M + sig_width * epsmid)
            sig_area = 0.22
            area = 10**(-3.42 + 0.90 * M + sig_area * epsmid)
        elif mech is Mechanism.R:
            sig_length = 0.16
            length = 10**(-2.42 + 0.58 * M + sig_length * epsmid)
            sig_width = 0.15
            width = 10**(-1.61 + 0.41 * M + sig_width * epsmid)
            sig_area = 0.26
            area = 10**(-3.99 + 0.98 * M + sig_area * epsmid)
        elif mech is Mechanism.N:
            sig_length = 0.17
            length = 10**(-1.88 + 0.50 * M + sig_length * epsmid)
            sig_width = 0.12
            width = 10**(-1.14 + 0.35 * M + sig_width * epsmid)
            sig_area = 0.22
            area = 10**(-2.78 + 0.82 * M + sig_area * epsmid)
        elif mech is Mechanism.A:
            sig_length = 0.16
            length = 10**(-2.44 + 0.59 * M + sig_length * epsmid)
            sig_width = 0.15
            width = 10**(-1.01 + 0.32 * M + sig_width * epsmid)
            sig_area = 0.24
            area = 10**(-3.49 + 0.91 * M + sig_area * epsmid)
        else:
            raise TypeError('Unsupported value of \'mech\'')
    elif rup_dim_model is MagScaling.S14:
        # Somerville (2014) model:
        #     - No length or width
        #     - No mechanism dependence
        sig_area = 0.3
        area = 10**(M - 4.25 + sig_area * epsmid)
        length = None
        sig_length = None
        width = None
        sig_width = None
    elif rup_dim_model == MagScaling.HB08:
        # Hanks and Bakun (2008)
        # These are the equations reported in the paper:
        #     M =       log10(A) + 3.98   for A <= 537 km^2 w/ se=0.03
        #     M = 4/3 * log10(A) + 3.07   for A >  537 km^2 w/ se=0.04
        # Using them is not so straight-forward beacuse we need to compute
        # the area from magnitude. Of course, this gives a different result
        # than if the equations were regressed for A as a function of M,
        # although since the equations were derived theoretically, this may
        # not be so bad.
        #
        # The inverted equation is simple enough:
        #     log10(A) =      M - 3.98    for M <= 6.71
        #     log10(A) = 3/4*(M - 3.07))  for M > 6.71
        #
        # The standard deviations are a little trickier.
        # First, convert standard errors of M to standard deviations of M:
        # (by my count, n=62 for A<=537, and n=28 for A>537)
        #     0.03*sqrt(62) = 0.236       for M <= 6.71
        #     0.04*sqrt(28) = 0.212       for M > 6.71
        # And convert to standard deviations of log(A) using the partial
        # derivatives
        #     dM/d(log10(A)) = 1          for M <= 6.71
        #     dM/d(log10(A)) = 3/4        for M >  6.71
        # So
        #     0.236*1   = 0.236 (pretty close to WC94)
        #     0.212*3/4 = 0.159 (seems a bit low...)
        if M > 6.71:
            sig_area = 0.236
            area = 10**(3 / 4 * (M - 3.07) + sig_area * epsmid)
        else:
            sig_area = 0.159
            area = 10**((M - 3.98) + sig_area * epsmid)
        length = None
        sig_length = None
        width = None
        sig_width = None
    elif rup_dim_model == MagScaling.Sea10_interface:
        # Strasser et al. (2010), coefficients for interface events
        sig_length = 0.18
        length = 10**(-2.477 + 0.585 * M + sig_length * epsmid)
        sig_width = 0.173
        width = 10**(-0.882 + 0.351 * M + sig_width * epsmid)
        sig_area = 0.304
        area = 10**(-3.49 + 0.952 * M + sig_area * epsmid)
    elif rup_dim_model == MagScaling.Sea10_slab:
        # Strasser et al. (2010), coefficients for slab events
        sig_length = 0.146
        length = 10**(-2.35 + 0.562 * M + sig_length * epsmid)
        sig_width = 0.067
        width = 10**(-1.058 + 0.356 * M + sig_width * epsmid)
        sig_area = 0.184
        area = 10**(-3.225 + 0.89 * M + sig_area * epsmid)
    else:
        raise TypeError('Unsupported value of \'rup_dim_model\'')
    return length, sig_length, width, sig_width, area, sig_area


def compute_epsilon(neps, trunc):
    """
    Compute midpoints and probabilities of epsilon bins.

    Args:
        neps (int): The number of steps to integrate from -trunc to +trunc.
            Larger numbers increase the accuracy of the result, but take
            longer to run.
        trunc (float): For the integration in area (or length and width), trunc
            is the truncation of the normal distribution (in units of sigma).
    Returns:
        tuple: epsilon midpoints, their probabilities, bin width.
    """
    # Need to assume a truncation level for normal distribution
    eps = np.linspace(-trunc, trunc, neps + 1)
    epsmid = 0.5 * (eps[1:] + eps[:-1])
    peps = (norm.cdf(eps[1:]) - norm.cdf(eps[:-1]))

    # define delta epsilons to normalize probabilities
    d_eps = 2 * trunc / neps
    epsfac = np.trapz(peps, dx=d_eps)

    if neps > 1:
        peps /= epsfac
    return epsmid, peps, d_eps
