# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
from os import path

cdir = path.abspath(path.dirname(__file__))

setup(
      name='snowmicrotoolset',
      version='0.0.1',
      description=u'Python tools for processing SnowMicroPen\u00AE observations',
      long_description=u'Python tools for processing SnowMicroPen\u00AE observations. For more information on the SnowMicroPen\u00AE, see https://www.slf.ch/en/services-and-products/research-instruments/snowmicropen-r-smp4-version.html',
      url='https://github.com/m9brady/SMP_to_CSV',
      author='Climate Processes Section',
      author_email='',
      license='Open Government Licence  Canada',
      classifiers=[
              'Development Status :: 2 - Pre-Alpha',
              'Intended Audience :: Science/Research',
              'Topic :: Scientific/Engineering',
              'License :: Other/Proprietary License',
              'Programming Language :: Python :: 2',
              'Programming Language :: Python :: 2.7',
              ],
      keywords='snow microstructure smp snowmicropen snowmicropentrometer',
      packages=find_packages(where=cdir, exclude=['contrib','docs','scratch','.git*','.spyproject']),
      # run-time dependencies
      install_requires=[
              'numpy',
              'scipy',
              'pandas',
              'matplotlib',
              'geojson',
              'folium'
              ],
        python_requires='~=2.7',	
        # additional dependencies
        extras_require={
                'dev': ['check-manifest'],
                'test': ['coverage'],
                }
    )