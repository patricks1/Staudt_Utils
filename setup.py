#!/usr/bin/env python
from distutils.core import setup
from staudt_utils import __version__

setup(name = 'staudt_utils',
      version = __version__,
      description = 'Miscellaneous utilities that I find useful',
      author='Patrick Staudt',
      author_email='patrickstaudt1@gmail.com',
      url='',
      platforms=['*nix'],
      license='GPL',
      requires = ['numpy','astropy','IPython','math'],
      provides = ['staudt_utils'],
      packages = ['staudt_utils']
      )
