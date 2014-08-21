calcSlopeDegrees
================

A Python script to calculate slope from a DEM.

The horizontal spacing of the DEM can be in metres or degrees. If in degrees the pixel size in m will be calculated for each pixel based on latitude.

If the DEM is noisy an option is available to fit a plane over a window of pixels and calculate the slope from this.

There are three implementations of the slope calculation, the fastest is chosen based on available libraries.

* Pure Python
* Numba JIT compiled Python (non-plane fitting version only)
* Fortran (wrapped using f2py)

This script was developed within the Microwave Systems, Sensors and Imaging Lab (MiXIL) at the University of Southern California (USC) as part work supported through NASA's Making Earth System Data Records for Use in Research Environments (MEaSUREs) Program. 

Building
---------

The script requires the following libraries:
* NumPy
* RIOS
* Numba (optional - if Fortran module is not available)
* ATLAS (optional - required on Linux to build the Fortran module)

To build the Fortran module and install the script use:

```
python setup.py build
python setup.py install --prefix=/install/dir
```

Usage
-----
DEM horizontal spacing in degrees:

```
calcSlopeDegrees.py --spacing_degrees inDEM.kea outSlope.kea
```

Fit plane over 9 x 9 window:
```
calcSlopeDegrees.py --plane_ls --window_size 9 inDEM.kea outSlope.kea
```

