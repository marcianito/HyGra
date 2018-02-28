# HyGra: hydro-gravitational modeling with R
setups, calculations, pre- and postprocessing

## Description

This R-package is a collection of function which serve as a basis for hydro-gravitational projects.
It evolved and grew throughout my PhD and hopefully helps others when dealing with this kind of data and projects.
It provides a structure for the following aspects, concerning hydro-gravitational projects and calculations:

* model domain preparation using DEMs (absolute UTM coordinates or relative coordinates)
* generating gravity component grids (circle, rectangular)
* remove effects of gravimeter pillar (circle, rectangular)
* determine conversion factor between gravity <> water (based on topography)
* forward modeling of gravity from water storage changes

Furthermore, this package consts the basis for further research projects which include coding in R, namely:

* [Estimating and correcting gravity data for the umbrella effect](https://github.com/marcianito/UmbrellaEffect)
* [Hydro-gravitational infiltration modeling](https://github.com/marcianito/gravityInf)
* [Synthetic uncertainty analysis of standard gravity corrections](https://github.com/marcianito/gravitySynth)

I am interested and open for further development and colaboration.
For bug fixes, comments or further development please contact: mreich@posteo.de or use github issues.

## Installation

1. Start R
2. Install package via devtools: 
`devtools::install_github("marcianito/HyGra")`

3. load package: 
`library(HyGra)`

## Dependencies

* r-base version 3.3.1
* following R-packages: devtools, dplyr, ggplot2, gstat, ptinpoly, raster, reshape2, viridis, zoo
* system libraries for rgdal and devtools

in debian install using: 
`apt-get update && apt-get install libgdal-dev libproj-dev libssl-dev`

