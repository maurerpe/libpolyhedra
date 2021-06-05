# libpolyhedra
C library to analyze and manipulate polyhedra with triangular faces.  Tested on linux.  Should work on any POSIX compliant OS.  Windows support untested.

## Features
* Reading and writing `.obj` and binary `.stl` files.
* Calculation of volume, center of mass, and inertia tensor
* Simplification
* Convex Hull
* Cut with a plane
* Convex decomposition

## Linux Build Instructions
```
$ autoreconf -i
$ ./configure
$ make
$ sudo make install
```
See `INSTALL` for more information.

## Algorithms
Includes algorithms from various publications.  See `PAPERS` for a list.  The authors of those papers were **not** involved with libpolyhedra and do not necessarily endorse it.
