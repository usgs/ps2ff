[![Build Status](https://travis-ci.org/usgs/ps2ff.svg?branch=master)](https://travis-ci.org/usgs/ps2ff)
[![codecov](https://codecov.io/gh/usgs/ps2ff/branch/master/graph/badge.svg)](https://codecov.io/gh/usgs/ps2ff)


# ps2ff
Produce approximated finite fault distances and variance corrections given point source 
information (e.g., Repi (epcentral distance) to Rjb (Joyner-Boore distance) or Rrup 
(closest distance to rupture).


DISCLAIMER:
This software has been approved for release by the U.S. Geological Survey (USGS).
Although the software has been subjected to rigorous review, the USGS reserves the
right to update the software as needed pursuant to further analysis and review. No
warranty, expressed or implied, is made by the USGS or the U.S. Government as to
the functionality of the software and related material nor shall the fact of release
constitute any such warranty. Furthermore, the software is released on condition
that neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from its authorized or unauthorized use.


Running the Programs
--------------------

The programs are Python (either 2.7 or 3.5) and require Numpy and Scipy. The easiest
way to get Numpy and Scipy is to install Anaconda. The programs will run on multiple
processors (cores) on machines so equipped.

To run the programs, do:

    % program.py config_file.ini

Program Details
---------------

*RjbMeanVar.py* 

Produces tables of Repi-to-Rjb ratios and variances. Example config
file "test_Rjb.ini". 

The parameters are:

 - **NP** - The number of processors (cores) to use. Minimum 1.

 - **filebase** - The base name of the output file (see output file naming convention below).

 - **datadir** - The directory into which the output files are written.

 - **rup_dim_model** - If 'WC94' use Wells and Coppersmith (1994) model, otherwise use Somerville
(2014).

 - **mech** - The rupture mechanism. One of "A" (all mechanisms), "SS" (strike-slip), "N" (normal), 
or "R" (reverse).

 - **LW** - (True or False) Use separate rupture Length and Width distributions rather than a 
rupture area distribution with a fixed aspect ratio.

 - **AR** - The Aspect Ratio (Length/Width) of the rupture. The aspect ratio will be maintained
until the rupture width spans the seismogenic zone, after which only the rupture 
length will increase.

 - **min_seis_depth** - The minimum seismogenic depth (km, integer).

 - **max_seis_depth* - The maximum seismogenic depth (km, integer).

 - **mindip_deg** - The minimum rupture dip in degrees (0 min, 90 max, integer).

 - **maxdip_deg** - The maximum rupture dip in degrees (0 min 90 max, integer).

 - **ndip** - The number of discrete steps in dip to use in the integration. Larger numbers
increase the accuracy of the result, but take longer to run.

 - **ntheta** - The number of discrete steps in theta to use in the integration. Larger numbers
increase the accuracy of the result, but take longer to run.

 - **nxny** - The number of discrete steps in x and y to use in the integration. Larger numbers
increase the accuracy of the result, but take longer to run. (Min 2)

 - **trunc** - For the integration in area (or length and width), trunc is the truncation
of the normal distribution (in units of sigma).

 - **neps** - The number of steps to integrate from -trunc to +trunc. Larger numbers
increase the accuracy of the result, but take longer to run.

 - **minmag** - The minimum magnitude for which to compute results (float).

 - **maxmag** - The maximum magnitude for which to compute results (float).

 - **dmag** - The size of the steps from minmag to maxmag (float).

 - **minepi** - The minimum epicentral distance for which to compute results (float).

 - **maxepi** - The maximum epicentral distance for which to compute results (float).

 - **nepi** - The number of steps from minepi to max epi. The steps will be uniformly sized
in log space.

*RrupMeanVar.py*

Produces tables of Repi-to-Rrup ratios and variances. Example config
file "test_Rrup.ini". 

The parameters are the same as for RjbMeanVar.py, above, with the
addition of:

 - **nz** - The number of steps in depth for Ztor. For any given rupture width and dip in the 
integration, Ztor ranges from (max_seis_depth - width * sin(dip)) to min_seis_depth.

*RrupRjbMeanVar_SingleEvent.py* 

Produces tables of Repi-to-Rrup and Repi-to-Rjb ratios and variances as a function of 
backazimuth for a particular earthquake magnitude and hypocentral depth. 
Example config file "test_single.ini". 

The parameters NP, datadir, rum_dim_model, mech, AR, ndip, mindip_deg, maxdip_deg, ntheta, 
nxny, neps, trunc, minepi, maxepi, nepi, min_seis_depth, and max_seis_depth, are the same 
as for *RjbMeanVar.py*, above, with the addition of:
 - **M** - The earthquake magnitude (float).
 - **zhyp** - The hypocentral depth of the earthquake (float).

Output File Naming Convention
------------------------------

The file names take the form:

    Rjb_S14_mechA_ar1p0_seis0_15_Ratios.csv

where:
 - "Rjb" is the the *filebase* parameter in the configuration file.
 - "S14" is the *rup_dim_model* (either "S14" or "WC94").
 - "mechA" specifies the rupture mechanism parameter *mech*, where "A" is one of "A", 
    "SS", "N", or "R".
 - "ar1p0" is the aspect ratio specified with the *AR* parameter, where the decimal point
    is replaced with the letter 'p'.
 - "seis0_15" is the range min/max seismogenic depths (in this case 0 to 15 km).
 - "Ratios" is either "Ratios" or "Var" specifying whether the file contains Repi to Rjb
    (or Rrup) ratios or variances.

Output File Format
------------------

Each output file starts with six header lines (each beginning with "#") specifying
the processing parameters. This is followed by a line (comma-separated) providing the
column headers. The first column, "Repi_km", is the epicentral distance. The following
columns R(magnitude) ("R" for "ratio") or V(magnitude) ("V" for "variance) provide the
values for a given Repi and magnitude. The table is intended for bi-variate interpolation,
linear in magnitude and logarithmic in distance.
The ratios are Rjb (or Rrup) to Repi. 

