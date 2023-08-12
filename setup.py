from setuptools import setup,find_packages

from glob import glob

import os
import sys
import platform


setup(name='modelham',
    version='0.1',
    description='A QM package for model Hamiltonian.',
    url='http://github.com/ZhaoYilin/modelham',
    author='Yilin Zhao',
    author_email='zhaoyilin1991@gmail.com',
    license='MIT',
    packages=find_packages())
