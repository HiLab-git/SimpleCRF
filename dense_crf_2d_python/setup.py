import sys
import numpy 
from distutils.core import setup
from distutils.extension import Extension

version = sys.version[0]
package_name = 'crf_max_flow'

module_name = 'dense_crf'

module = Extension(module_name,
                   include_dirs=[numpy.get_include(), './densecrf/include', './densecrf/external/liblbfgs/include'],
                   sources=['wrap_py{0:}.cpp'.format(version),
                            'dense_crf.cpp',
                            'dense_crf_core.cpp',
                            './densecrf/src/densecrf.cpp',
                            './densecrf/src/labelcompatibility.cpp',
                            './densecrf/src/objective.cpp',
                            './densecrf/src/optimization.cpp',
                            './densecrf/src/pairwise.cpp',
                            './densecrf/src/permutohedral.cpp',
                            './densecrf/src/unary.cpp',
                            './densecrf/src/util.cpp',
                            './densecrf/external/liblbfgs/lib/lbfgs.c'
                            ])

setup(name=package_name,
      ext_modules = [module])


# to build, run python stup.py build or python setup.py build_ext --inplace
# to install, run python setup.py install