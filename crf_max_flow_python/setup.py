import sys
import numpy
from distutils.core import setup
from distutils.extension import Extension

version = sys.version[0]
package_name = 'crf_max_flow'
module_name1 = 'max_flow'
module_name2 = 'max_flow3d'

wrap_source_2d = "wrap_2d_py{0:}.cpp".format(version)
wrap_source_3d = "wrap_3d_py{0:}.cpp".format(version)
module1 = Extension(module_name1,
                    include_dirs = [numpy.get_include()],
                    sources = ['max_flow.cpp', 'maxflow-v3.0/graph.cpp', 'maxflow-v3.0/maxflow.cpp', wrap_source_2d])
module2 = Extension(module_name2,
                    include_dirs = [numpy.get_include()],
                    sources = ['max_flow3d.cpp', 'maxflow-v3.0/graph.cpp', 'maxflow-v3.0/maxflow.cpp', wrap_source_3d])
setup(name=package_name,
      ext_modules = [module1, module2])


# to build, run python stup.py build or python setup.py build_ext --inplace
# to install, run python setup.py install
