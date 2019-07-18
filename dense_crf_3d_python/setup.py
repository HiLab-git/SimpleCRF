import sys
import numpy
from distutils.core import setup
from distutils.extension import Extension

version = sys.version[0]
package_name = 'crf_max_flow'
module_name  = 'dense_crf3d'

module = Extension(module_name,
                    include_dirs = [numpy.get_include(), 'densecrf_3d/include', 'densecrf_3d/external/liblbfgs/include'],
                    sources = ['wrap_py{0:}.cpp'.format(version),
                               './dense_crf3d.cpp', 
                               './dense_crf3d_util.cpp', 
                               './densecrf_3d/src/densecrf.cpp', 
                               './densecrf_3d/src/labelcompatibility.cpp',
                               './densecrf_3d/src/objective.cpp', 
                               './densecrf_3d/src/optimization.cpp',
                               './densecrf_3d/src/pairwise.cpp',
                               './densecrf_3d/src/permutohedral.cpp',
                               './densecrf_3d/src/unary.cpp',
                               './densecrf_3d/src/util.cpp',
                               './densecrf_3d/external/liblbfgs/lib/lbfgs.c'])
setup(name=package_name,
      ext_modules = [module])


# to build, run python stup.py build or python setup.py build_ext --inplace
# to install, run python setup.py install
