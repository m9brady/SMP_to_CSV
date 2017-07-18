from setuptools import setup, find_packages
from codecs import open
from os import path

cdir = path.abspath(path.dirname(__file__))

setup(
	name='snowmicrotoolset',
	version='0.0.1',
	description='Python tools for processing SMP observations',
	long_description='',
	url='https://github.com/m9brady/SMP_to_CSV',
	author='Climate Processes Section',
	author_email='',
	license='Open Government Licence – Canada',
	
	classifiers=[
		'Development Status :: 2 - Pre-Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering',
		'License :: Other/Proprietary License',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
	],
	
	keywords='snow smp snowmicropentrometer',
	
	packages=find_packages(exclude=['contrib','docs','tests']),
	
	# run-time dependencies
	install_requires=[
		'numpy',
		'scipy'
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