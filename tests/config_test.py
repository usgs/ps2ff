#!/usr/bin/env python

import os.path
import tempfile
import textwrap
import sys
import pytest
import copy

from configobj import ConfigObj

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))

from validate import ValidateError

import ps2ff.config as config
from ps2ff.constants import DistType, MagScaling, Mechanism

def test_config():

    #
    # Break the config
    #
    ctest = {}
    ctest['min_seis_depth'] = 10.
    ctest['max_seis_depth'] = 5.
    with pytest.raises(ValidateError):
        config.check_config(ctest)
    ctest['max_seis_depth'] = 20.

    ctest['mindip_deg'] = 10.
    ctest['maxdip_deg'] = 5.
    with pytest.raises(ValidateError):
        config.check_config(ctest)

    ctest['maxdip_deg'] = 10.
    ctest['ndip'] = 5
    with pytest.raises(ValidateError):
        config.check_config(ctest)
    ctest['ndip'] = 1

    ctest['minmag'] = 10.
    ctest['maxmag'] = 5.
    with pytest.raises(ValidateError):
        config.check_config(ctest)
    ctest['maxmag'] = 15.

    ctest['minepi'] = 10.
    ctest['maxepi'] = 5.
    with pytest.raises(ValidateError):
        config.check_config(ctest)

    ctest['maxepi'] = 10.
    ctest['nepi'] = 5.
    with pytest.raises(ValidateError):
        config.check_config(ctest)

    ctest['nepi'] = 1.

    config.check_config(ctest)

    #
    # magScalingType()
    #
    with pytest.raises(ValidateError):
        res = config.magScalingType(['1', '2'])
    with pytest.raises(ValidateError):
        res = config.magScalingType('NotAThing')
    res = config.magScalingType('WC94')
    assert res is MagScaling.WC94

    #
    # mechType()
    #
    with pytest.raises(ValidateError):
        res = config.mechType(['1', '2'])
    with pytest.raises(ValidateError):
        res = config.mechType('NotAThing')
    res = config.mechType('R')
    assert res is Mechanism.R



if __name__ == '__main__':
    test_config()
