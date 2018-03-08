# PSR Geometry

----
## Unreleased

* Generalised the Bstep() function, for arbitrary vector fields.
* Include retardation in the calculation of the emission point.

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

