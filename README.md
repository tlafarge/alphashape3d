# alphashape3d

The package __alphashape3d__ presents the implementation in R of the __alpha-shape__ of a finite set of points in the __three-dimensional__ space. This geometric structure generalizes the convex hull and allows to recover the shape of non-convex and even non-connected sets in 3D, given a random sample of points taken into it. Besides the computation of the alpha-shape, the package __alphashape3d__ provides users with functions to compute the volume of the alpha-shape, identify the connected components and facilitate the three-dimensional graphical visualization of the estimated set. 

To install:

* the latest released version from CRAN: `install.packages("alphashape3d")`
* the latest development version without vignette:                              `install_github("alphashape3d","tlafarge",build_vignettes = FALSE)` (require  "devtools")
* the latest development version: `install_github("alphashape3d","tlafarge")` (require "devtools")

Authors: 

* Beatriz Pateiro-Lopez <beatriz.pateiro@usc.es>
* Thomas Lafarge <thomas.lafarge@nist.gov>

License:

* GPL-2
