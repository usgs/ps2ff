#!/usr/bin/env python3

import numpy as np
from scipy.stats import norm


def dimensions_from_magnitude(M, rup_dim_model, neps, trunc, mech='A'):
    """
    Compute dimensions of rupture from magnitude for a specified
    magnitude scaling relation.

    Args:
        M (float): Magnitude.
        rup_dim_model (str): String indicating the model for compputing the
            rupture dimensions from magnitude. Supported values are:
                - 'WC94'
                - 'S14'
                - 'HB08'
        neps (int): The number of steps to integrate from -trunc to +trunc.
            Larger numbers increase the accuracy of the result, but take
            longer to run.
        trunc (float): For the integration in area (or length and width), trunc
            is the truncation of the normal distribution (in units of sigma).
        mech (str): Optional string indicating earthquake mechanism, used by
            some of the models. Anything other than 'R', 'N', 'SS', or 'A'
            (the default).

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

    if rup_dim_model == 'WC94':
        # Use mech to get either M-A or (M-W) and (M-R) from Wells and
        # Coppersmith.
        if mech == "SS":
            sig_length = 0.15
            length = 10**(-2.57 + 0.62 * M + sig_length * epsmid)
            sig_width = 0.14
            width = 10**(-0.76 + 0.27 * M + sig_width * epsmid)
            sig_area = 0.22
            area = np.power(10, -3.42 + 0.90*M +
                            sig_area * epsmid)
        if mech == "R":
            sig_length = 0.16
            length = 10**(-2.42 + 0.58 * M + sig_length * epsmid)
            sig_width = 0.15
            width = 10**(-1.61 + 0.41 * M + sig_width * epsmid)
            sig_area = 0.26
            area = np.power(10, -3.99 + 0.98*M +
                            sig_area * epsmid)
        if mech == "N":
            sig_length = 0.17
            length = 10**(-1.88 + 0.50 * M) + sig_length * epsmid
            sig_width = 0.12
            width = 10**(-1.14 + 0.35 * M + sig_width * epsmid)
            sig_area = 0.22
            area = np.power(10, -2.78 + 0.82*M +
                            sig_area * epsmid)
        if mech == "A":
            sig_length = 0.16
            length = 10**(-2.44 + 0.59 * M + sig_length * epsmid)
            sig_width = 0.15
            width = 10**(-1.01 + 0.32 * M + sig_width * epsmid)
            sig_area = 0.24
            area = np.power(10, -3.49 + 0.91*M +
                            sig_area * epsmid)
    elif rup_dim_model == 'S14':
        # Somerville (2014) model:
        #     - No length or width
        #     - No mechanism dependence
        sig_area = 0.3
        area = np.power(10, M - 4.25 + sig_area * epsmid)
        length = None
        sig_length = None
        width = None
        sig_width = None
    else:
        raise Exception('Unsupported value of \'rup_dim_model\'')
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
    eps = np.linspace(-trunc, trunc, neps+1)
    epsmid = 0.5*(eps[1:] + eps[:-1])
    peps = (norm.cdf(eps[1:]) - norm.cdf(eps[:-1]))

    # define delta epsilons to normalize probabilities
    d_eps = 2 * trunc / neps
    epsfac = np.trapz(peps, dx=d_eps)

    if neps > 1:
        peps /= epsfac
    return epsmid, peps, d_eps
