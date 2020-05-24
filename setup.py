# -*- coding: utf-8 -*-
from os import path
from re import search

from setuptools import find_packages, setup

cdir = path.abspath(path.dirname(__file__))

with open(path.join(cdir, "snowmicrotoolset", "__init__.py")) as f:
    version = search(r'__version__\s*=\s*"(\S+)"', f.read()).group(1)

setup(
    name='snowmicrotoolset',
    version=version,
    description=u'Python tools for processing SnowMicroPen\u00AE observations',
    long_description=u'Python tools for processing SnowMicroPen\u00AE observations. ' +\
        u'For more information on the SnowMicroPen\u00AE instrument, see ' +\
        u'https://www.slf.ch/en/services-and-products/research-instruments/snowmicropen-r-smp4-version.html',
    url='https://github.com/m9brady/SMP_to_CSV',
    author='Climate Processes Section',
    author_email='',
    license='Open Government Licence - Canada',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: Other/Proprietary License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='snow microstructure smp snowmicropen snowmicropentrometer',
    packages=find_packages(
        where=cdir, 
        exclude=['contrib','docs','scratch','.git*','.spyproject']
    ),
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'tqdm',
        'click'
    ],
    python_requires='>2.6, !=3.0, !=3.1, !=3.2, !=3.3, !=3.4, !=3.5, <3.9',
    # additional dependencies
    extras_require={
        ':python_version == "2.7"': ['futures'],
        'dev': [
            'ipython',
        ],
    },
    entry_points='''
    [console_scripts]
    smp_to_csv=snowmicrotoolset.scripts.cli:cli
    '''
)
