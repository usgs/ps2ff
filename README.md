ps2ff
=====
Produce approximated finite fault distances and variance corrections given
point source information (e.g., Repi (epcentral distance) to Rjb (Joyner-Boore
distance) or Rrup (closest distance to rupture).

[![Build Status](https://travis-ci.org/usgs/ps2ff.svg?branch=master)](https://travis-ci.org/usgs/ps2ff)
[![codecov](https://codecov.io/gh/usgs/ps2ff/branch/master/graph/badge.svg)](https://codecov.io/gh/usgs/ps2ff)


Prerequisites and Installation
------------------------------
These programs run in either Pything 2.7 or 3.5, and require a number of
packages. The dependency list is given in `install.sh`. The easiest way to
install this code is to run `install.sh` in OSX or Linux. It will install
miniconda (if a version of `conda` is not already installed) and then all
the dependencies in a virtual environment named `ps2ff`. To use this
environment after installation, type
```
source activate ps2ff
```

Running the Programs
--------------------
The primary program is `run_ps2ff`, which takes only one argument
```
ps2ff config_file.ini
```
There are example config files in `<repository>/config`. 

Output Tables
-------------
The 'tables' directory contains example results for some generic seismological
assumptions. The output file name convension is easiest to describe with an
example:
```
Rjb_S14_mechA_ar1p0_seis0_15_Ratios.csv
```
where:
 - "Rjb" is the the `what` parameter in the configuration file.
 - "S14" is the selected `rup_dim_model`.
 - "mechA" specifies the rupture mechanism parameter `mech`, where "A" can
   be one of "A", "SS", "N", or "R".
 - "ar1p0" is the aspect ratio specified with the `AR` parameter, where the
   decimal point is replaced with the letter 'p'.
 - "seis0_15" is the range min/max seismogenic depths (in this case 0 to 15
   km).
 - "Ratios" is either "Ratios" or "Var" specifying whether the file contains
   Repi to Rjb (or Rrup) ratios or variances.

Each output table starts with six header lines (each beginning with "#")
specifying the processing parameters. This is followed by a line
(comma-separated) providing the column headers. The first column, "Repi_km",
is the epicentral distance. The following columns R(magnitude) ("R" for
"ratio") or V(magnitude) ("V" for "variance) provide the values for a given
Repi and magnitude. The table is intended for bi-variate interpolation, linear
in magnitude and logarithmic in distance. The ratios are Rjb (or Rrup) to Repi.


Program Details
---------------

`run_ps2ff` produces tables of Repi-to-Rjb ratios and variances. Example config
file "test_Rjb.ini". The parameters in the config file are:

 - `NP` The number of processors (cores) to use. Minimum 1.

 - `filebase` The base name of the output file (see output file naming
   convention below).

 - `datadir` - The directory into which the output files are written. If
   unspecified, it uses "./data".

 - `rup_dim_model` String to select the magnitude scaling relationship.
   Currently supported values are:
        - 'WC94' for Wells and Coppersmith (1994),
        - 'S14' for Somerville (2014).

 - `mech` The rupture mechanism, only used by some scaling relationships:
        - 'A' for all/unknown mechanisms)
        - 'SS' for strike-slip,
	- 'N' for normal,
	- 'R' for reverse.

 - `LW` Boolean for whether to separately select rupture length and width
   distributions, otherwise select the rupture area and compute length and
   width from it and an assumed aspect ratio. 

 - `AR` Aspect ratio (Length/Width) of the rupture. The aspect ratio is
   maintained until the rupture width spans the seismogenic zone, after
   which only the rupture length will increase.

 - `min_seis_depth` The minimum seismogenic depth (km).

 - `max_seis_depth` The maximum seismogenic depth (km).

 - `mindip_deg` The minimum rupture dip in degrees (0 min, 90 max).

 - `maxdip_deg` The maximum rupture dip in degrees (0 min 90 max).

 - `ndip` The number of integration steps in dip.

 - `ntheta` The number of integration steps in theta.

 - `nxny` The number of integration steps in x and y (minimum is 2).

 - `trunc` For the integration in area (or length and width), this is the 
   truncation of the normal distribution (in standard deviation units).

 - `neps` The number of integration steps for area (or length and width)
   from -trunc to +trunc. Larger numbers increase the accuracy of the result,
   but take longer to run.

 - `minmag` The minimum magnitude for which to compute results.

 - `maxmag` The maximum magnitude for which to compute results.

 - `dmag` The size of the steps from minmag to maxmag.

 - `minepi` The minimum epicentral distance for which to compute results.

 - `maxepi` - The maximum epicentral distance for which to compute results.

 - `nepi` The number of steps from minepi to max epi. The steps will be 
    uniformly sized in log space.

 - `nz` The number of integration steps in depth for Ztor. For any given
   rupture width and dip in the integration, Ztor ranges from 
   `(max_seis_depth - width * sin(dip))` to `min_seis_depth`. Only used for
   if `what='Rrup'`. 


`RrupRjbMeanVar_SingleEvent.py` roduces tables of Repi-to-Rrup and Repi-to-Rjb
ratios and variances as a function of backazimuth for a particular earthquake
magnitude and hypocentral depth. Example config file "test_single.ini".

The parameters are the same as for `run_ps2ff`, with the addition of:
 - `M` The earthquake magnitude.
 - `zhyp` The hypocentral depth of the earthquake.
