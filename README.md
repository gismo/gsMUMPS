# gsMUMPS

## Installation
1. Obtain the MUMPS files from [http://mumps-solver.org](http://mumps-solver.org)
2. Inside the `/path/to/MUMPS.X.Y.Z` directory, see the `INSTALL` file for instructions.
	- **EXAMPLE**: For common Linux systems, this line is enough:
```
$ cp Make.inc/Makefile.debian.PAR ./Makefile.inc
```
	- *NOTE*: Make sure to build the shared library of MUMPS, done by:
```
$ make allshared
```
3. Inside `gismo`, turn on this module:
```
cmake . -DGISMO_OPTIONAL="gsMUMPS;<other modules>"
```

