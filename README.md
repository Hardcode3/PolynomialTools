# PolynomialTools
## Description

A simple C++ / CMake project to deal with n'th order polynomial resolution and other simple polynomial operations.

The project is currently under development and is delivered with unit tests.

The main feature is the resolution of polynomials of degree n by using the method of the companion matrix applied to its monic form. The eigen values of this matrix correspond to its roots and both real and imaginary roots are found by using this method.
Depending on the needs, you can pick real, imaginary of every roots of the polynomial.
This resolution can be useful in many topics of applied mathematics.

Results have been compared to numpy.roots() function from matplotlib (python library), giving the exact same output.  

## Dependencies

The dependency is limited to the library [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)

## Installation

In order to use the library:
- Embed the needed files in your project
- Or clone the whole project as a git submodule in your existing project

Eigen can be:
- An existing target in your host repository
- An external target (vcpkg...)
- A local target created by cloning the repository when setting the flag POLYNOMIALS_CLONE_SUBMODULE_EIGEN to ON

You are free to chose the intergration in your project.

FindPackage is not suported on this library.

## Sources
The algorithm has been writen using:
- The [Wikipedia page](https://en.wikipedia.org/wiki/Companion_matrix)
- [R. A. Horn & C. R. Johnson, Matrix Analysis. Cambridge, UK: Cambridge University Press, 1999, pp. 146-7.](https://anandinstitute.org/pdf/Roger_A.Horn.%20_Matrix_Analysis_2nd_edition(BookSee.org).pdf)
