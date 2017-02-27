### __**HyGra** is an R-package for hydro-gravitational modeling__ ###

This package was created during my PhD in hydrogravimetry.
It contains a catalog of functions which I collected and developed during the work on the anaylsis for my publications.
Its focus is on combining gravity and hydrological modeling, data analysis, pre- and postprocessing as also plotting.
This package is certainly not complete but unique in the R-universe in combining this type of data analysis.
With the functions and workflows included, it is feasible to carry out a hydro-gravimetric study (terrestrial gravimetry) from the first setup, through data measurements in the field,
setting up a hydrological model and combining data again for joint analysis and figure generation (plots).
Many things will have to be adjusted to individual needs.
I am interested and open for further development and colaboration.
Enjoy!

#### __Installation:__ ####

    * load package "devtools":

          library(devtools)

    * install package via

          install_bitbucket("marciano/HyGra")
          
    * load package 
    
          library(hygra)

#### __Documentation:__ ####

Documentation can be accessed via the Vignete or using 

    ?hygra

#### __Main functionality:__ ####

    * read, organized and filter various field measurement data
    * gravity grid preparation
    * gravity signal calculations (1D, 2D & 3D)
    * HYDRUS pre- and postprocessing
    * some hydrological parameter estimation functions
    * cross-correlations of different hydrological time series
    * hydro-gravimetric infiltration model (simulating a sprinkling experiment)

### Contact ###

For bugs and specific problems please open an issue.

Contact via email: mreich@posteo.de
