from enum import Enum


class AutoName(Enum):
    def __init__(this):
        this._value_ = this.name


"""
These enum members hold their name (as a string) as their value.
For example: DistType.Rjb.value == 'Rjb'
"""


class DistType(AutoName):
    """
    Type of distance measure to compute.
       - 'Rjb' is for the Joyner-Boore distance.
       - 'Rrup' is for closest distance to rupture.
    """
    Rjb = ()
    Rrup = ()


class MagScaling(AutoName):
    """
    Magnitude scaling relationship.
       - 'WC94' is for Wells and Coppersmith (1999)
       - 'HB08' is for Hanks and Bakun (2008)
       - 'S14' is for Somerville et al. (2014)
       - 'Sea10_interface' is for Strasser et al. (2010), for interface events
       - 'Sea10_slab' is for Strasser et al. (2010), for intraplate events
    """
    WC94 = ()
    HB08 = ()
    S14 = ()
    Sea10_interface = ()
    Sea10_slab = ()


class Mechanism(AutoName):
    """
    Source mechanism.
       - 'A' for all (mechanism is unknown/unspecified)
       - 'R' for reverse
       - 'N' for normal
       - 'SS' for strike-slip
    """
    A = ()
    R = ()
    N = ()
    SS = ()
