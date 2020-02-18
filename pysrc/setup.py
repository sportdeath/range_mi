from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import os

# Add install locations
include_dirs = ["/usr/local/include", np.get_include()]
library_dirs = ["/usr/local/lib"]

extensions = [
        Extension(
            "grid_mi",
            ["grid_mi.pyx"],
            libraries=["range_mi"],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            language='c++',
            extra_compile_args=["-std=c++11", "-Ofast"],
            extra_link_args=["-std=c++11"]
            )]

setup(
    ext_modules=cythonize(extensions)
)
