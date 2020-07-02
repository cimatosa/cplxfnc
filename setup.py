from setuptools import setup
from setuptools.extension import Extension
from ctypes.util import find_library
from Cython.Build import cythonize

import pathlib
import subprocess

author      = "Richard Hartmann"
authors     = [author]
description = 'special functions with complex arguments - wrapped for python'
name        = 'cplxfnc'
version     = '0.1.1'

def try_to_install(names):
    if (names[0] == 'arb' or names[0] == 'flint'):
        return install_arb()

    print("no install procedure found for {}".format(names))
    return False

def install_arb():
    r = input("Do you want to download, build and install arb 2.15.0 from https://github.com/fredrik-johansson/arb? [y/n]")
    if r != 'y':
        return False

    dir = pathlib.Path('.').absolute()

    print("install dependencies libgmb and libmpfr")
    print("---------------------------------------")
    subprocess.check_call('sudo apt-get install libgmp-dev libmpfr-dev', shell=True)

    print("install flint2")
    print("--------------")

    print("flint2 is installed from source, this requires to remove libflint-* from the package manager")
    r = input("do you want to continue [y/n]")
    if r != 'y':
        return False

    subprocess.check_call('sudo apt-get remove libflint-*', shell=True)

    subprocess.check_call('wget -O flint2.zip https://github.com/fredrik-johansson/flint2/archive/trunk.zip', shell=True)
    subprocess.check_call('unzip flint2.zip', shell=True)
    subprocess.check_call('cd flint2-trunk && ./configure && make && sudo make install', shell=True)

    print("install arb")
    print("-----------")

    fname = 'arb-2.15.tar.gz'

    subprocess.check_call('wget -O {} https://github.com/fredrik-johansson/arb/archive/2.15.0.tar.gz'.format(fname), shell=True)
    subprocess.check_call('tar -zxvf {}'.format(fname), shell=True)
    subprocess.check_call('cd arb-2.15.0 && ./configure && make && sudo make install', shell=True)

    subprocess.check_call('sudo /sbin/ldconfig', shell=True)

    return True



def check_for_lib(*names):
    lib_found = None
    for n in names:
        l = find_library(n)
        if l is not None:
            lib_found = n
            print("found c library {}".format(l))
            break

    if lib_found is None:
        print("shared library '{}' not found".format(names))
        if not try_to_install(names):
            raise RuntimeError("Error!")


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