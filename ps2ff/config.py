import os.path
import pkg_resources

from ps2ff.constants import DistType, MagScaling, Mechanism

from configobj import ConfigObj, flatten_errors
from validate import Validator, ValidateError

def get_configspec():
    """
    Returns the full path to the configspec.conf file.

    Returns:
        (str): The path to the configspec.

    """
    return os.path.join(pkg_resources.resource_filename('ps2ff', 'data'),
                        'configspec.conf')

def get_custom_validator():
    """
    Returns a validator suitable for validating the ps2ff config file.

    Returns:
        (:class:`Validator`): A Validator object.

    """
    fdict = {
        'distType': distType,
        'magScalingType': magScalingType,
        'mechType': mechType,
    }
    validator = Validator(fdict)
    return validator

def config_error(config, results):
    """
    Parse the results of a ConfigObj validation and print the errors.
    Throws a RuntimeError exception  upon completion if any errors or 
    missing sections are encountered.

    Args:
        config (ConfigObj): The ConfigObj instance representing the 
            parsed config.
        results (dict): The dictionary returned by the validation of
        the 'config' arguments.

    Returns:
        (Nothing): Nothing

    """
    errs = 0
    for (section_list, key, _) in flatten_errors(config, results):
        if key is not None:
            print('The "%s" key in the section "%s" failed validation' %
                    (key, ', '.join(section_list)))
            errs += 1
        else:
            print('The following section was missing:%s ' %
                    ', '.join(section_list))
            errs += 1
    if errs:
        raise RuntimeError('There %s %d %s in configuration.' %
                    ('was' if errs == 1 else 'were', errs,
                     'error' if errs == 1 else 'errors'))

def check_config(config):
    """
    Checks the config for sanity.
    """
    if config['min_seis_depth'] >= config['max_seis_depth']:
        print('Config error: min_seis_depth must be < max_seis_depth')
        raise ValidateError()
    if config['mindip_deg'] > config['maxdip_deg']:
        print('Config error: mindip_deg must be <= maxdip_deg')
        raise ValidateError()
    if config['mindip_deg'] == config['maxdip_deg'] and config['ndip'] != 1:
        print('Config error: mindip_deg == maxdip_deg, but ndip != 1')
        raise ValidateError()
    if config['minmag'] > config['maxmag']:
        print('Config error: minmag must be <= maxmag')
        raise ValidateError()
    if config['minepi'] > config['maxepi']:
        print('Config error: minepi must be <= maxepi')
        raise ValidateError()
    if config['minepi'] == config['maxepi'] and config['nepi'] != 1:
        print('Config error: minepi == maxepi, but nepi != 1')
        raise ValidateError()

def distType(value):
    """
    Checks that the value is in the DistType enum and returns 
    the appropriate enum.
    """
    if not isinstance(value, str):
        print('Config error: Unknown value "%s" for "what"' % (value))
        raise ValidateError()
    for thing in DistType:
        if value == thing.value:
            return thing
    print('Config error: Invalid value "%s" for "what"' % (value))
    print('Config error: "what" must be in %s' % 
         ([x.value for x in list(DistType)]))
    raise ValidateError()

def magScalingType(value):
    """
    Checks that the value is in the magScalingType enum and returns 
    the appropriate enum.
    """
    if not isinstance(value, str):
        print('Config error: Unknown value "%s" for "rup_dim_model"' % (value))
        raise ValidateError()
    for thing in MagScaling:
        if value == thing.value:
            return thing
    print('Config error: Invalid value "%s" for "rup_dim_model"' % (value))
    print('Config error: "rup_dim_model" must be in %s' % 
         ([x.value for x in list(MagScaling)]))
    raise ValidateError()

def mechType(value):
    """
    Checks that the value is in the mechType enum and returns 
    the appropriate enum.
    """
    if not isinstance(value, str):
        print('Config error: Unknown value "%s" for "mech"' % (value))
        raise ValidateError()
    for thing in Mechanism:
        if value == thing.value:
            return thing
    print('Config error: Invalid value "%s" for "mech"' % (value))
    print('Config error: "mech" must be in %s' % 
         ([x.value for x in list(Mechanism)]))
    raise ValidateError()
