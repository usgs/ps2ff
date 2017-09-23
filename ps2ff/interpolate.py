"""
PS2FF class for converting point distances (epicentral or hypocentral)
to equivalent average finite rupture distances (Rjb or Rrup). Based
upon:
    Thomposn, E. M., C. B. Worden (2018). Estimating rupture distances
    without a rupture Bull. Seism. Soc. Am. (in press).
"""
import os.path
import re
import glob
import pkg_resources

import numpy as np
import pandas as pd
import scipy.interpolate as spint

from ps2ff.constants import DistType, MagScaling, Mechanism


class PS2FF(object):
    @classmethod
    def fromParams(cls, dist_type=DistType.Rjb, mag_scaling=MagScaling.WC94,
                   mechanism=Mechanism.A, AR=1.7, min_seis_depth=0,
                   max_seis_depth=20):
        """
        Create a PS2FF object from a set of parameters. Parameters must
        combine into a file name of an existing table in the ps2ff
        tables. The file name will take the form:

        <dist_type>_<mag_scaling>_mech<mechanism>_ar<AR>_seis<min_seis_depth>_<max_seis_depth>_Ratios.csv

        where the decimal point in the aspect ratio will be replaced
        with the letter 'p'.

        Args:
            dist_type (DistType): One of the DistType enum members.
                Typically DistType.Rjb (default) or DistType.Rrup. See
                ps2ff.constants for a complete list.
            mag_scaling (MagScaling): One of the MagScaling enum
                members.See ps2ff.constants for a complete list.
                The default is MagScaling.WC94.
            mechanism (Mechanism): A mechanism from the Mechanism
                enum. See ps2ff.constants for a complete list. The
                default is Mechanism.A.
            AR (float): The aspect ratio for the rupture computations.
                Typical values are 1.7 (default) and 1.0; tables may
                not exist for other values.
            min_seis_depth (int): The depth (km) to the top of the
                seismogenic zone. This is typically 0 (default);
                tables for other values may not exist.
            max_seis_depth (int): The depth (km) to the bottom of the
                seismogenic zone. Typical values are 15 (for ACR
                regions) and 20 (default, for SCR regions).

        Returns:
            (PS2FF): An object of the PS2FF class initialized with
            the tables corresponding to the selected parameters.
        """  # noqa
        filebase = '%s_%s_mech%s_ar%.1f_seis%d_%d' % (dist_type.value,
                                                      mag_scaling.value,
                                                      mechanism.value, AR,
                                                      int(min_seis_depth),
                                                      int(max_seis_depth))
        filebase = filebase.replace('.', 'p')
        rfile = filebase + '_Ratios.csv'

        datadir = pkg_resources.resource_filename('ps2ff', 'tables')
        filepath = os.path.join(datadir, rfile)
        if not os.path.isfile(filepath):
            print('Error: could not find file with the supplied parameters')
            print('File not found: %s' % (rfile))
            print('The available files are:\n%s' %
                  '\n'.join(cls.getValidFiles()))
            raise ValueError('Unknown file %s' % (rfile))
        return cls.fromFile(rfile)

    @classmethod
    def fromFile(cls, rfile):
        """
        Create a PS2FF object from a file specification. The file must
        exist in the ps2ff tables.

        Args:
            rfile (str): A file name (base file name only, not a full
                  path) corresponding to one of the tables in the
                  ps2ff.data resource. The file should be the
                  "Ratios" file; the "Var" file name is derived
                  from the Ratios file.

        Returns:
            (PS2FF): An object of the PS2FF class initialized with
            the tables corresponding to the rfile argument.
        """
        this = cls()
        datadir = pkg_resources.resource_filename('ps2ff', 'tables')
        filepath = os.path.join(datadir, rfile)
        if not os.path.isfile(filepath):
            print('File not found: %s' % (rfile))
            print('The available files are:\n%s' %
                  '\n'.join(cls.getValidFiles()))
            raise ValueError('Unknown file %s' % (rfile))

        r2r_ratios_tbl = pd.read_csv(filepath, comment='#')
        r2r_cols = r2r_ratios_tbl.columns[1:]
        mag_list = []
        for column in (r2r_cols):
            if re.search('R\d+\.*\d*', column):
                magnitude = float(re.findall('R(\d+\.*\d*)', column)[0])
                mag_list.append(magnitude)
        mag_list = np.array(mag_list)

        r2r_dist_name = r2r_ratios_tbl.columns[0]
        dist_list = np.log(np.array(r2r_ratios_tbl[r2r_dist_name]))

        r2r_ratios_grid = r2r_ratios_tbl.values[:, 1:]

        this.r2r_ratios_obj = spint.RectBivariateSpline(dist_list, mag_list,
                                                        r2r_ratios_grid,
                                                        kx=1, ky=1)

        varfile = filepath.replace('_Ratios', '_Var')
        r2r_var_tbl = pd.read_csv(varfile, comment='#')

        r2r_var_grid = r2r_var_tbl.values[:, 1:]
        this.r2r_var_obj = spint.RectBivariateSpline(dist_list, mag_list,
                                                     r2r_var_grid,
                                                     kx=1, ky=1)
        this.rfile = rfile
        this.vfile = os.path.basename(varfile)
        return this

    def r2r(this, r, M):
        """
        Convert point distances to the equivalent average finite rupture
        distances, based on the parameters specified when creating this
        object.

        Args:
            r (numpy.ndarray): An array of point distances (typically
                epicentral distance) in km.
            M (numpy.ndarray): An array (the same shape as r) of
                magnitudes.

        Returns:
            (numpy.ndarray): An array the same shape as r, with
                distances converted to average finite rupture distance.
        """
        return r * this.r2r_ratios_obj.ev(np.log(r), M)

    def rat(this, r, M):
        """
        Return ratios needed to convert point distances to the equivalent
        average finite rupture distances, based on the parameters specified
        when creating this object.

        Args:
            r (numpy.ndarray): An array of point distances (typically
                epicentral distance) in km.
            M (numpy.ndarray): An array (the same shape as r) of
                magnitudes.

        Returns:
            (numpy.ndarray): An array the same shape as r, with
                the ratios (multipliers) needed to convert r to
                average finite rupture distance.
        """
        return this.r2r_ratios_obj.ev(np.log(r), M)

    def var(this, r, M):
        """
        Return the additional variance from the uncertainty in
        point distances vis a vis finite rupture distance.

        Args:
            r (numpy.ndarray): An array of point distances (typically
                epicentral distance) in km.
            M (numpy.ndarray): An array (the same shape as r) of
                magnitudes.

        Returns:
            (numpy.ndarray): An array the same shape as r, containing
            the additional variance from the uncertainty in finite
            rupture distance for the distances in r.
        """
        return this.r2r_var_obj.ev(np.log(r), M)

    def files(this):
        """
        Returns the table files that were used to construct this
        object.

        Args:
            None

        Returns:
            (str, str): A tuple of the ratio file and variance file.
        """
        return this.rfile, this.vfile

    @classmethod
    def getValidFiles(cls):
        """
        Get a list of the valid conversion tables in the ps2ff
        package.

        Args:
            None

        Returns:
            A list of file names corresponding to the available tables
            in the ps2ff package.
        """
        datadir = pkg_resources.resource_filename('ps2ff', 'tables')
        validfiles = glob.glob(os.path.join(datadir, '*_Ratios.csv'))
        return [os.path.basename(x) for x in validfiles]
