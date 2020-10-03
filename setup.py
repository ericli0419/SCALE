#!/usr/bin/env python
"""
# Author: Yuzhe Li
# Created Time : Fri Sep 18 14:21:41 CST 2020

# File Name: setup.py
# Description:

"""

from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='scale-scanpy',
      version='0.0.2',
      description='Single-Cell ATAC-seq Analysis via Latent feature Extraciton',
      packages=find_packages(),

      author='Yuzhe Li',
      author_email='liyuzhezju@gmail.com',
      url='https://github.com/ericli0419/SCALE.git',
      scripts=['SCALE.py'],
      install_requires=requirements,
      python_requires='>3.6.0',

      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.7',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: POSIX :: Linux',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
     ],
     )
