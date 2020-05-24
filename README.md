# SMP_to_CSV
Converts SMP measurement files (```*.pnt```) into CSV containing all data points (Depth and Penetration Force) plus a multi-line header containing time, latitude, longitude and other site-specific information. A depth/force profile quicklook is also generated as a PNG file.

## Disclaimer
The SMP class methods ```retrieve_header```,  ```extract_data``` and general function ```xcorr``` are from the excellent [SnowMicroPyn](https://sourceforge.net/projects/pyntreader) tool by Sascha Grimm. SnowMicroPyn is licensed under the [GNU General Public License version 3.0 (GPLv3)](https://sourceforge.net/directory/license:gplv3/).

The ```detect_peaks``` and associated ```_plot``` helper functions are from the [Biomechanics and Motor Control](https://github.com/demotu/BMC/blob/master/functions/detect_peaks.py) repo by Marcos Duarte. BMC is licensed under the [MIT Licence](https://github.com/demotu/BMC/blob/master/LICENSE.txt)

## Installation
*The following assumes you have a working Python 2.7 or 3.6/3.7/3.8 environment which includes `pip`*


### Option A: Direct from GitHub
Use the install-from-github functionality from `pip` to install from the `master` branch:
```
pip install git+https://github.com/m9brady/SMP_to_CSV.git@master
```

### Option B: Manual Clone
Make a local clone of this repo, navigate to where it is located and run the following command to add the module to your current python environment:
```
pip install .
```

## Package Usage

### Option A: Batch processing SMP to CSV/PNG

When `snowmicrotoolset` is installed, a new utility is added to your Python scripts path: `smp_to_csv`. Simply run this tool with no arguments to get some help documentation on what it needs to function:

```
$ smp_to_csv

***************************
** No input folder given **
***************************

Usage: smp_to_csv [OPTIONS]

  Convert SLF .PNT files to CSV/PNG formats

Options:
  -i, --input_folder TEXT   Location of input .smp files
  -o, --output_folder TEXT  Location of output csv/png files  [default: ./outdata]
  --version                 Show the version and exit.
  -h, --help                Show this message and exit.
```

### Option B: Interactively processing PNT files

You should now be able to import the SMP class in a python interpreter or script using the following:
```
>>> from snowmicrotoolset import SMP
>>> p = SMP('path/to/pnt_file.pnt')
```

## License
This dataset is licensed under the [Open Government License of Canada](http://open.canada.ca/en/open-government-licence-canada)
and is subject to the [Copyright Act of Canada](http://laws-lois.justice.gc.ca/eng/acts/C-42/index.html). Additional information can be found at the [Government of Canada's Open Government portal](http://open.canada.ca)
