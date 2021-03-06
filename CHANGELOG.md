# PSR Geometry

## Unreleased

* Fit alpha and zeta for given main pulse and interpulse PPAs.

----
## v1.2.4 (2018-03-23)

* Fixed bug in the calculation of the polarisation position angle.

----
## v1.2.3 (2018-03-22)

* Fixed bug in the calculation of the A field.
* Added capability of calc\_field\_test to make sure A was perpendicular to V.

----
## v1.2.2 (2018-03-20)

* Corrected erroneous Vstep() function, and renamed it traj\_step().

----
## v1.2.1 (2018-03-19)

* Added a Vstep() function, for stepping along velocity field lines.
* Implemented psr\_trajectory, for plotting particle trajectories.

----
## v1.2.0 (2018-03-18)

* Implemented a "fitwidth" function (and program) for matching emission heights to pulse widths.

----
## v1.1.12 (2018-03-18)

* Restricted emission point range to <90% of light cylinder radius.
* Debugged support for the -i option (opposite beam) in psr\_polangle and psr\_emit.

----
## v1.1.11 (2018-03-17)

* Include retardation in the calculation of the emission point.

----
## v1.1.10 (2018-03-17)

* Fixed a convergence bug in the convergence of the emission point.

----
## v1.1.9 (2018-03-16)

* Fixed important bug in psr\_polangle (still in beta)

----
## v1.1.8 (2018-03-13)

* Add program psr\_polangle for plotting PA curves (beta)

----
## v1.1.7 (2018-03-13)

* Added support for Powell optimisation as an alternative to Nelder-Mead optimisation.
* Added support for the 'elevator' method (recommended)

----
## v1.1.6 (2018-03-11)

* Added program psr\_cost\_function for evaluating the cost functions at a grid of x,y,z coordinates.

----
## v1.1.5 (2018-03-11)

* Added functions for calculating the line of sight, and converting between position angles and beam opening angles.
* Added a function for finding the approximate emission point, using a simplified magnetic field geometry.
* Half-fixed bugs relating to finding the emission point in the full Arendt geometry.

## v1.1.4 (2018-03-09)

* Added an output format option "BLF" to psr\_fields.

----
## v1.1.3 (2018-03-08)

* Added a function to find the emission point on a last open field line corresponding to a given rotation phase.
* Added a function, norm\_dot(), for calculating the normalised dot product of two vectors (i.e points).
* Added the ability to plot a custom range of z's in psr\_field

----
## v1.1.2 (2018-03-07)

* Added a program for plotting magnetic field lines
* Added an option to footpoint() and farpoint() functions to normalise the output coordinates to the light cylinder radius.
* Added an option to footpoint() and farpoint() functions to stop if the point gets too far out.

----
## v1.1.1 (2018-03-06)

* Added a function for finding a magnetic field lines most extreme point.
* Extended the functionality of "psr\_fields" to print out a whole lot of different quantities.

----
## v1.1.0 (2018-03-02)

* Added/updated man pages for all the functions in this library.
* Removed deprecated Bfield() and Vfield() functions (which were subsumed into calc\_fields())

----
## v1.0.4 (2018-03-02)

* Fixed known bugs with the calculation of the acceleration field

----
## v1.0.3 (2018-03-02)

* Added a small, utility program, "psr\_fields", which prints out the B and V fields as columns of numbers

----
## v1.0.2 (2018-03-02)

* Added a calc\_fields() function that does the job of Bfield() and Vfield() and also calculates the acceleration vector

----
## v1.0.1 (2018-02-21)

* Changed the "angle" struct and functions to "psr\_angle"

----
## v1.0.0 (2018-02-21)

* Initial set of features implemented

