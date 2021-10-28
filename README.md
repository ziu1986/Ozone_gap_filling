# Ozone gap filling
Reynolds decomposition for gap filling of ozone data.

## Installation
-------------
Download and install e.g. anaconda and run:

`$ conda env create -f environment.yml`

Activate environment

`$ conda activate ozone_gap_filling`

To use ipython once the environment is installed, the python module `decorator` has to be installed by hand.

`$ conda install -n ozone_gap_filling decorator`

## Release note
-------------
### v.1.0

This Repository contains the functions and workflow to use Reynolds decomposition on ozone monitoring data.
All functions were optimized for python 2.7.

## Data & Preprocessing
------------------------
Ozone monitoring data as obtained from ebas.nilu.no.
Copernicus regional ozone reanalysis as obtained from [Copernicus regional air monitoring](https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-europe-air-quality-forecasts?tab=form) and 
the annual mean was computed using cdo/nco tool kits.

## Workflows
Currently there is only the "reconstruction" workflow defined.

## Running
-----------
You can run directly from the command line

`$ python ozone_gap_filling.py`

or use the ipython interpreter

```ipython
$ ipython

In [1]: execfile('ozone_gap_filling.py')
```