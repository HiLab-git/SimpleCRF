from distutils.core import setup
from distutils.extension import Extension

package_name = 'dense_crf'
module_name1 = 'dense_crf'

module1 = Extension(module_name1,
                    include_dirs = ['/Users/guotaiwang/anaconda/lib/python2.7/site-packages/numpy/core/include','include', 'external/liblbfgs/include'],
                    sources = ['dense_crf.cpp', 'src/densecrf.cpp', 'src/labelcompatibility.cpp','src/objective.cpp', 'src/optimization.cpp','src/pairwise.cpp','src/permutohedral.cpp','src/unary.cpp','src/util.cpp','external/liblbfgs/lib/lbfgs.c'])
setup(name=package_name,
      ext_modules = [module1])


# to build, run python stup.py build or python setup.py build_ext --inplace
# to install, run python setup.py install
