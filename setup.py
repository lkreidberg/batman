from __future__ import print_function
import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import numpy as np
import os
import shutil
import tempfile


class CustomBuildExt(build_ext):
    def build_extensions(self):
        ce = self.compiler
        if ce.compiler_type == "unix":
            if self.compiler_has_openmp():
                for e in self.extensions:
                    e.extra_compile_args += ["-fopenmp", "-std=c99"]
                    e.libraries += ["gomp"]
            else:
                for e in self.extensions:
                    e.extra_compile_args += ["-std=c99"]
        super().build_extensions()

    def compiler_has_openmp(self):
        """Check if the compiler supports OpenMP"""
        test_code = "#include <omp.h>\nint main() { omp_get_num_threads(); return 0; }"
        return self.try_compile(test_code, extra_postargs=["-fopenmp"])

    def try_compile(self, code, extra_postargs=None):
        """Attempt to compile a test program"""
        tmp_dir = tempfile.mkdtemp(prefix="tmp-setuptools")
        file_name = os.path.join(tmp_dir, "test.c")
        with open(file_name, "w") as file:
            file.write(code)
        try:
            self.compiler.compile(
                [file_name], output_dir=tmp_dir, extra_postargs=extra_postargs
            )
            return True
        except setuptools.distutils.errors.CompileError:
            return False
        finally:
            shutil.rmtree(tmp_dir)


extensions = [
    Extension("batman._nonlinear_ld", ["c_src/_nonlinear_ld.c"]),
    Extension("batman._quadratic_ld", ["c_src/_quadratic_ld.c"]),
    Extension("batman._uniform_ld", ["c_src/_uniform_ld.c"]),
    Extension("batman._logarithmic_ld", ["c_src/_logarithmic_ld.c"]),
    Extension("batman._exponential_ld", ["c_src/_exponential_ld.c"]),
    Extension("batman._custom_ld", ["c_src/_custom_ld.c"]),
    Extension("batman._power2_ld", ["c_src/_power2_ld.c"]),
    Extension("batman._rsky", ["c_src/_rsky.c"]),
    Extension("batman._eclipse", ["c_src/_eclipse.c"]),
]

setup(
    name="batman-package",
    version="2.5.2",
    author="Laura Kreidberg",
    author_email="laura.kreidberg@gmail.com",
    url="https://github.com/lkreidberg/batman",
    packages=["batman"],
    license="GNU GPLv3",
    description="Fast transit light curve modeling",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
    ],
    include_dirs=[np.get_include()],
    install_requires=["numpy"],
    setup_requires=["wheel"],
    extras_require={
        "matplotlib": ["matplotlib"],
    },
    ext_modules=extensions,
    cmdclass={"build_ext": CustomBuildExt},
)
