from distutils.core import setup
from distutils.extension import Extension
import numpy as np

occultnl = Extension('batman.occultnl', ['c_src/occultnl.c'], extra_compile_args = ['-fopenmp'], extra_link_args=['-lgomp'])
occultquad = Extension('batman.occultquad', ['c_src/occultquad.c'], extra_compile_args=['-fopenmp'], extra_link_args=['-lgomp'])
occultuniform = Extension('batman.occultuniform', ['c_src/occultuniform.c'], extra_compile_args = ['-fopenmp'], extra_link_args=['-lgomp'])
rsky = Extension('batman.rsky', ['c_src/rsky.c'])

setup(	name='batman', 
	version='1.0', 
	packages =['batman'],
	include_dirs = [np.get_include()],
	ext_modules=[occultnl, occultquad, occultuniform, rsky]
)
