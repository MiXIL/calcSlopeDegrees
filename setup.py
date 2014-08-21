#!/usr/bin/env python
"""
Setup script for calc slope script.
"""

import platform
from numpy.distutils.core import setup, Extension

# Build with ATLAS on linux
if platform.system() == 'Linux':
    liblist = ['lapack','f77blas','cblas','latlas']
# Else try to build with LAPACK (will work on OS X)
else:
    liblist = ['lapack']

slope = Extension('slope',
                    include_dirs = ['/usr/include','/usr/local/include'],
                    libraries = ['lapack'],
                    library_dirs = ['/usr/lib','/usr/local/lib'],
                    sources = ['slope.f90'])

setup(
  name='slopeDegrees',
  version = '1.0',
  description = 'Script to calculate slope from a DEM with a horizontal spacing in degrees',
  author = 'Daniel Clewley',
  author_email = 'daniel.clewley@gmail.com',
  url = 'http://mixil.usc.edu/',
  ext_modules = [slope],
  scripts=['calcSlopeDegrees.py']
)
