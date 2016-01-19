from __future__ import print_function
from distutils.core import setup
from distutils.extension import Extension
import numpy as np
from distutils.ccompiler import new_compiler
import os
import sys
import tempfile

"""
Check for OpenMP based on
https://github.com/MDAnalysis/mdanalysis/tree/develop/package/setup.py
retrieved 06/15/15
"""
def detect_openmp():
	"""Does this compiler support OpenMP parallelization?"""
	compiler = new_compiler()
	print("Checking for OpenMP support... ")
	hasopenmp = hasfunction(compiler, 'omp_get_num_threads()')
	needs_gomp = hasopenmp
	if not hasopenmp:
		compiler.add_library('gomp')
	hasopenmp = hasfunction(compiler, 'omp_get_num_threads()')
	needs_gomp = hasopenmp
	if hasopenmp: print("Compiler supports OpenMP")
	else: print( "Did not detect OpenMP support.")
	return hasopenmp, needs_gomp

def hasfunction(cc, funcname, include=None, extra_postargs=None):
	# From http://stackoverflow.com/questions/
	#            7018879/disabling-output-when-compiling-with-distutils
	tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
	devnull = oldstderr = None
	try:
		try:
			fname = os.path.join(tmpdir, 'funcname.c')
			f = open(fname, 'w')
			if include is not None:
				f.write('#include %s\n' % include)
			f.write('int main(void) {\n')
			f.write('    %s;\n' % funcname)
			f.write('}\n')
			f.close()
			# Redirect stderr to /dev/null to hide any error messages
			# from the compiler.
			# This will have to be changed if we ever have to check
			# for a function on Windows.
			devnull = open('/dev/null', 'w')
			oldstderr = os.dup(sys.stderr.fileno())
			os.dup2(devnull.fileno(), sys.stderr.fileno())
			objects = cc.compile([fname], output_dir=tmpdir, extra_postargs=extra_postargs)
			cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
		except Exception as e:
			return False
		return True
	finally:
		if oldstderr is not None:
			os.dup2(oldstderr, sys.stderr.fileno())
		if devnull is not None:
			devnull.close()

#checks whether OpenMP is supported
has_openmp, needs_gomp = detect_openmp()
parallel_args = ['-fopenmp'] if has_openmp else []
parallel_libraries = ['gomp'] if needs_gomp else []

_nonlinear_ld = Extension('batman._nonlinear_ld', ['c_src/_nonlinear_ld.c'], extra_compile_args = parallel_args, libraries = parallel_libraries) 
_quadratic_ld = Extension('batman._quadratic_ld', ['c_src/_quadratic_ld.c'], extra_compile_args = parallel_args, libraries = parallel_libraries) 
_uniform_ld   = Extension('batman._uniform_ld', ['c_src/_uniform_ld.c'], extra_compile_args = parallel_args, libraries = parallel_libraries) 
_logarithmic_ld   = Extension('batman._logarithmic_ld', ['c_src/_logarithmic_ld.c'], extra_compile_args = parallel_args, libraries = parallel_libraries) 
_exponential_ld   = Extension('batman._exponential_ld', ['c_src/_exponential_ld.c'], extra_compile_args = parallel_args, libraries = parallel_libraries) 
_custom_ld   = Extension('batman._custom_ld', ['c_src/_custom_ld.c', 'c_src/_custom_intensity.c'], extra_compile_args = parallel_args, libraries = parallel_libraries) 
_rsky = Extension('batman._rsky', ['c_src/_rsky.c'])
_eclipse = Extension('batman._eclipse', ['c_src/_eclipse.c'])

setup(	name='batman-package', 
	version="2.2.0",
	author='Laura Kreidberg',
	author_email = 'laura.kreidberg@uchicago.edu',
	url = 'https://github.com/lkreidberg/batman',
	packages =['batman'],
	license = ['GNU GPLv3'],
	description ='Fast transit light curve modeling',
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering',
		'Programming Language :: Python'
		],
	include_dirs = [np.get_include()],
	install_requires = ['numpy'],
	ext_modules=[_nonlinear_ld, _quadratic_ld, _uniform_ld, _logarithmic_ld, _exponential_ld, _custom_ld, _rsky, _eclipse]
)
