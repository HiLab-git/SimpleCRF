from distutils.core import setup
from distutils.extension import Extension

package_name = 'crf_max_flow'
#module1 = Extension('dense_crf',
#                   include_dirs=['/Users/guotaiwang/anaconda/lib/python2.7/site-packages/numpy/core/include', './densecrf/include'],
#                   sources=['dense_crf.cpp', 'dense_crf_wrapper.cpp'],
#                   libraries = ['densecrf','lbfgs', 'objective', 'optimization','pairwise','permutohedral','unary','util'],
#                   library_dirs = ['./libs'])
module1 = Extension('dense_crf',
                   include_dirs=['/Users/guotaiwang/anaconda/lib/python2.7/site-packages/numpy/core/include', './densecrf/include','./densecrf/external/liblbfgs/include'],
                   sources=['dense_crf.cpp',
                            'dense_crf_wrapper.cpp',
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
      ext_modules = [module1])


# to build, run python stup.py build or python setup.py build_ext --inplace
# to install, run python setup.py install