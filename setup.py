from setuptools import setup
from setuptools.extension import Extension
from ctypes.util import find_library
from Cython.Build import cythonize

author      = "Richard Hartmann"
authors     = [author]
description = 'special functions with complex arguments - wrapped for python'
name        = 'cplxfnc'
version     = '0.1.1'

def check_for_lib(*names):
    lib_found = None
    for n in names:
        l = find_library(n)
        if l is not None:
            lib_found = n
            print("found c library {}".format(l))
            break

    if lib_found is None:
        raise RuntimeError("shared library '{}' not found".format(n))
    return lib_found

libraries = []
libraries.append(check_for_lib('m'))
libraries.append(check_for_lib('flint'))
libraries.append(check_for_lib('arb', 'flint-arb'))

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
                                       libraries = libraries,
                                       extra_compile_args=['-std=c++11'],
                                       language="c++"))
)