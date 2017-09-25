ps2ff
=====
[![Build Status](https://travis-ci.org/usgs/ps2ff.svg?branch=master)](https://travis-ci.org/usgs/ps2ff)
[![codecov](https://codecov.io/gh/usgs/ps2ff/branch/master/graph/badge.svg)](https://codecov.io/gh/usgs/ps2ff)

<img align="left" height="70" src="ps2ff/data/ps2ff.png">
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
 - "Rjb" is the the `what` parameter in the configuration file.
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
variances. The parameters in the config file are:

- `NP` The number of processors (cores) to use. Minimum 1.

- `datadir` The directory into which the output files are written. If
  unspecified, it uses `./data`.

- `rup_dim_model` String to select the magnitude scaling relationship.
    Currently supported values are:
  - WC94: Wells, D. L., & Coppersmith, K. J. (1994). New empirical
    relationships among magnitude, rupture length, rupture width, rupture area,
    and surface displacement, *Bulletin of the Seismological Society of
    America*, 84(4), 974-1002.
  - S14: Somerville, P. (2014). Scaling Relations between Seismic Moment and
    Rupture Area of Earthquakes in Stable Continental Regions, *PEER Report*
    2014/14.
  - HB08: Hanks, T. C. and Bakun, W. H. (2008). M-logA observations for
    recent large earthquakes, *Bulletin of the Seismological Society of
    America*, 98(1), 490-494.
  - Sea10_interface: Interface coefficients of Strasser, F. O., Arango,
    M. C., & Bommer, J. J. (2010). Scaling of the source dimensions of
    interface and intraslab subduction-zone earthquakes with moment magnitude,
    *Seismological Research Letters*, 81(6), 941-950.
  - Sea10_slab: Slab coefficients from the paper in previous bullet.

- `mech` The rupture mechanism, only used by some scaling relationships:

  - A: all/unknown mechanisms,
  - SS: strike-slip,
  - N: normal,
  - R: reverse.

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

- `maxepi` The maximum epicentral distance for which to compute results.

- `nepi` The number of steps from minepi to max epi. The steps will be
   uniformly sized in log space.

- `nz` The number of integration steps in depth for Ztor. For any given
  rupture width and dip in the integration, Ztor ranges from
  `(max_seis_depth - width * sin(dip))` to `min_seis_depth`. Only used for
  if `what='Rrup'`.


`run_ps2ff_single_event` produces tables of Rrup-to-Repi and Rjb-to-Repi
ratios and variances for a single event. This means that the magnitdue and
hypocentral depth are available, simplifying the integration. It optionally
tabulates the adjustment factors as a function of backazimuth. An example
configuration file for this program is given in
`tests/config/test_single.ini`.

The parameters are the same as for `run_ps2ff`, with the addition of:
- `M` The earthquake magnitude.
- `zhyp` The hypocentral depth of the earthquake.
- `bytheta` Tabulate factors for bins of theta.

