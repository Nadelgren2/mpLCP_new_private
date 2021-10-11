from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = "matrix_manipulation",
  cmdclass = {"build_ext": build_ext},
  ext_modules =
  [
    Extension("matrix_manipulation",
              ["matrix_manipulation.pyx"],
              extra_compile_args = ["-O0", "-fopenmp"],
              extra_link_args=['-fopenmp']
              )
  ]
)
