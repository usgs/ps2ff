
import numpy as np
import time as time

from ps2ff.integration_loops import single_event_inner_loop
from ps2ff.constants import MagScaling, Mechanism


# Use max_workers threading in model?

def single_event_adjustment(
        magnitude,
        hyp_depth,
        ar=1.7,
        mechanism=Mechanism.A,
        mag_scaling=MagScaling.WC94,
        n_repi=13,
        min_repi=0.1,
        max_repi=1000,
        nxny=3,
        n_theta=9,
        n_dip=3,
        min_dip=0,
        max_dip=np.pi/2.0,
        n_eps=3,
        trunc=2):
    """
    Method for getting event-specific point source distance
    adjustment factors. This does not integrate across magnitude
    or depth and so those values must be provided.

    Args:
        magnitude (float):
            Earthquake magnitude.
        hyp_depth (float):
            Hypocentral depth (km).
        ar (float):
            Aspect ratio (L/W).
        mechanism (Mechanism):
            A ps2ff.constants Mechanism instance.
        mag_scaling (MagScaling):
            A ps2ff.constants MagScaling instance.
        n_repi (int):
            Number of log-spaced Repi points to compute conversions.
        min_repi (float):
            Minimum Repi to compute conversions (km).
        max_repi (float):
            Maximum Repi to compute conversions (km).
        nxny (int):
            Number of integration steps in the x/y direction for floating,
            the rupture plane around the hypocenter.
        n_theta (int):
            Number of integration steps for theta. Default value of 9 is
            an angular step size of 45 deg since it goes from 0 to 360 deg.
        n_dip (int)
            Number of integration steps for dip. Default value of 3 gives
            a step size of 45 deg for the default range of 0 to 90 deg.
        min_dip (float):
            Minimum dip for integration (rad).
        max_dip (float):
            Maximum dip for integration (rad).
        n_eps (int):
            Number of integration steps for mag-area relationship.
        trunc (float):
            Truncation level for normal distribution of log area conditioned
            on magnitude.

    Retunrs:
        tuple: All arrays of length n_repi:

            - Repi (km).
            - Average Rjb conditioned on Repi, M, and Zhyp.
            - Average Rrup conditioned on Repi, M, and Zhyp.
            - Variance of Rjb conditioned on Repi, M, and Zhyp.
            - Variance of Rrup conditioned on Repi, M, and Zhyp.

    """
    # Check that mag_scaling and mechanism are the
    # correct type
    if not isinstance(mag_scaling, MagScaling):
        raise TypeError(
            'mag_scaling must bee a MagScaling instance')
    if not isinstance(mechanism, Mechanism):
        raise TypeError(
            'mechanism must bee a Mechanism instance')

    # Create conf dictionary.
    conf = {}
    conf['M'] = magnitude
    conf['zhyp'] = hyp_depth
    conf['AR'] = ar
    conf['mech'] = mechanism
    conf['rup_dim_model'] = mag_scaling

    conf['nxny'] = nxny

    conf['ndip'] = n_dip
    conf['mindip'] = min_dip
    conf['maxdip'] = max_dip

    conf['neps'] = n_eps
    conf['trunc'] = trunc

    conf['min_seis_depth'] = 0.0
    conf['max_seis_depth'] = 35.0  # not used?
    conf['bytheta'] = False

    repi = np.logspace(
        np.log10(min_repi), np.log10(max_repi), n_repi)

    Rrup_var = np.zeros_like(repi)
    Rrup_avg = np.zeros_like(repi)
    Rjb_var = np.zeros_like(repi)
    Rjb_avg = np.zeros_like(repi)

    for i in range(n_repi):
        Rrup_var[i], Rrup_avg[i], Rjb_var[i], Rjb_avg[i] = \
            single_event_inner_loop(conf, repi[i], ntheta=n_theta)

    return repi, Rjb_avg, Rrup_avg, Rjb_var, Rrup_var
