ps2ff
=====
[![Build Status](https://travis-ci.org/usgs/ps2ff.svg?branch=master)](https://travis-ci.org/usgs/ps2ff)
[![codecov](https://codecov.io/gh/usgs/ps2ff/branch/master/graph/badge.svg)](https://codecov.io/gh/usgs/ps2ff)

<img align="left" height="70" src="doc_source/_static/ps2ff_wide.png">
Produce approximated finite fault distances and variance corrections given
point source information (e.g., Repi (epcentral distance) to Rjb (Joyner-Boore
distance) or Rrup (closest distance to rupture).

<br><br>

Using the results (the API)
---------------------------

The command line programs (described below) can be used to generate new
distance adjustments. This package also includes a set of correction factors
for some common conditions (e.g., typical active crustal regions). These
can most easily be used with the `interpolate` module that contains the `PS2FF`
class, which enables the use of the tables for arbitrary magnitudes and
epicentral distance values. See the `ps2ff.interpolate` section of this
package's [documentation](https://usgs.github.io/ps2ff/).

Prerequisites and Installation
------------------------------
These programs run in Python 3.5 or later, and require a number of
packages, including the [conda](https://conda.io/docs/) package manager.
The dependency list is given in `install.sh`. The easiest way to install
this code is to run `install.sh` in OSX or Linux. It will install
miniconda (if a version of `conda` is not already installed) and then all
the dependencies in a virtual environment named `ps2ff`. To use this
environment after installation, type
```
source activate ps2ff
```

Running the Programs
--------------------
The primary program is `run_ps2ff`, which must be handed a configuraiton file
```
ps2ff config_file.ini
```
There are example configuration files in the `ps2ff/config` directory.

Output Tables
-------------
The `ps2ff/tables` directory contains example results for some generic seismological
assumptions. The output file name convension is easiest to describe with an
example:
```
Rjb_S14_mechA_ar1p0_seis0_15_Ratios.csv
```
where:
 - "Rjb" is the the `what` [command line argument](https://usgs.github.io/ps2ff/run_ps2ff.html).
 - "S14" is the selected `rup_dim_model`.
 - "mechA" specifies the rupture mechanism parameter `mech`, where "A" can
   be one of "A", "SS", "N", or "R".
 - "ar1p0" is the aspect ratio specified with the `AR` parameter, where the
   decimal point is replaced with the letter 'p'.
 - "seis0_15" is the range min/max seismogenic depths (in this case 0 to 15
   km).
 - "Ratios" is either "Ratios" or "Var" specifying whether the file contains
   Rjb- or Rrup-to-Repi ratios, or variances.

Each output table starts with six header lines (each beginning with `#`)
specifying the processing parameters. This is followed by a line
(comma-separated) providing the column headers. The first column, "Repi_km",
is the epicentral distance. The following columns "R(magnitude)" ("R" for
"ratio") or "V(magnitude)" ("V" for "variance) provide the values for a given
Repi and magnitude. The table is intended for bi-variate interpolation, linear
in magnitude and logarithmic in distance. The ratios are Rjb (or Rrup) to Repi.


Program Details
---------------

`run_ps2ff` produces tables of Rjb-to-Repi or Rrup-to-Repi ratios and
variances. Example configuration files may be found in ps2ff/config.

`run_ps2ff_single_event` produces tables of Rrup-to-Repi and Rjb-to-Repi
ratios and variances for a single event. This means that the magnitdue and
hypocentral depth are available, simplifying the integration. It optionally
tabulates the adjustment factors as a function of backazimuth. An example
configuration file for this program is given in
`tests/config/test_single.ini`.

See this package's documentation on the 
[command line interface](https://usgs.github.io/ps2ff/)
for more on usage and configuration options.
