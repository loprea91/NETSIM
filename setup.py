#!/usr/bin/env python

from distutils.core import setup
import glob

data_files = glob.glob('./src/**/*.*', recursive=True)

setup(name='pyNetsim',
      version='1.0',
      description='Python interface for NETSIM',
      author='',
      author_email='',
      url='',
      packages=['pyNetsim'],
      data_files=[('src', data_files)],
     )