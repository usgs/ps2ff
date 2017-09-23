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
    Rjb -> Joyner-Boore distance
    Rrup -> closest distance to rupture
    """
    Rjb = ()
    Rrup = ()


class MagScaling(AutoName):
    """
    Magnitude scaling relationship.
    WC94 -> Wells and Coppersmith, 1999
    HB08 -> Hanks and Bakun, 2008
    S14 -> Somerville et al., 2014
    Sea10_interface: Strasser et al. (2010), for interface events
    Sea10_slab: Strasser et al. (2010), for intraplate events
    """
    WC94 = ()
    HB08 = ()
    S14 = ()
    Sea10_interface = ()
    Sea10_slab = ()


class Mechanism(AutoName):
    """
    Source mechanism.
    A -> all (mechanism is unknown/unspecified)
    R -> reverse
    N -> normal
    SS -> strike-slip
    """
    A = ()
    R = ()
    N = ()
    SS = ()
