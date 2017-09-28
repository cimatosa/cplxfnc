from setuptools import setup
from setuptools.extension import Extension

from Cython.Build import cythonize

author      = "Richard Hartmann"
authors     = [author]
description = 'special functions with complex arguments - wrapped for python'
name        = 'cplxfnc'
version     = '0.1'

setup(
      name=name,
      author=author,
      url='https://github.com/cimatosa/cplxfnc',
      version=version,
      packages=[name],
      package_dir={name: name},
      license="MIT",
      description=description,
      long_description=description,
      ext_modules = cythonize(Extension("cplxfnc.cplxfnc_cyth",
                                       ["./cplxfnc/cplxfnc_cyth.pyx", "./cplxfnc_clib/cplxfnc.cpp"],
                                       libraries = ['m', 'flint', 'arb'],
                                       extra_compile_args=['-std=c++11'],
                                       language="c++"))
)