# gsMUMPS

## Installation
1. Obtain the MUMPS files from [http://mumps-solver.org](http://mumps-solver.org)
2. Inside the `/path/to/MUMPS.X.Y.Z` directory, see the `INSTALL` file for instructions.
	- **EXAMPLE**: For common Linux systems, this line is enough:
	```
	$ cp Make.inc/Makefile.debian.PAR ./Makefile.inc
	```
3. To avoid modifying the `LD_LIBRARY_PATH` upon every use, modify the following lines in `./Makefile.inc`
    ```
    # Adapt/uncomment RPATH_OPT to avoid modifying
    # LD_LIBRARY_PATH in case of shared libraries
    # RPATH_OPT = -Wl,-rpath,/path/to/MUMPS_x.y.z/lib/
    ```
    By uncommenting the `RPATH_OPT` and inserting the correct `/path/to/MUMPS.X.Y.Z`
    ```
    # Adapt/uncomment RPATH_OPT to avoid modifying
    # LD_LIBRARY_PATH in case of shared libraries
    RPATH_OPT = -Wl,-rpath,/path/to/MUMPS.X.Y.Z/lib/ # <------- HERE
    ```
4. After modifying the `Makefile.inc`, make the shared libraries of MUMPS
    ```
    # Clean previous builds (optional)
    $ make clean
    # Build
    $ make allshared
    ```
5. Inside `gismo`, turn on this module:
    ```
    $ cmake . -DGISMO_OPTIONAL="gsMUMPS;<other modules>" -DMUMPS_DIR=/path/to/MUMPS_X.Y.Z
    ```
6. Build `gismo` as usual
    ```
    $ make gismo
    ```
