# SMP_to_CSV
Converts SMP measurement files (```*.pnt```) into CSV containing all data points (Depth and Penetration Force) plus a multi-line header containing time, latitude, longitude and other site-specific information. A depth/force profile quicklook is also generated as a PNG file.

## Disclaimer
The SMP class methods ```retrieve_header```,  ```extract_data``` and general function ```xcorr``` are from the excellent [SnowMicroPyn](https://sourceforge.net/projects/pyntreader) tool by Sascha Grimm. SnowMicroPyn is licensed under the [GNU General Public License version 3.0 (GPLv3)](https://sourceforge.net/directory/license:gplv3/).

The ```detect_peaks``` and associated ```_plot``` helper functions are from the [Biomechanics and Motor Control](https://github.com/demotu/BMC/blob/master/functions/detect_peaks.py) repo by Marcos Duarte. BMC is licensed under the [MIT Licence](https://github.com/demotu/BMC/blob/master/LICENSE.txt)

## Requirements and Usage
*Assumes you have Python 2.7 installed with pip, numpy, scipy, pandas and matplotlib modules*

### Option A: Importing the SMP class definition

Make a local clone of this repo, navigate to where it is located and run the following command to add the module to your current python environment:
```
pip install .
```
You should now be able to import the SMP class in a python interpreter or script using the following:
```
from snowmicrotoolset import SMP
```

### Option B: Running the tool standalone

1. Create two folders ```indata``` and ```outdata``` in the same directory as ```snowmicrotoolset/__init__.py```

2. Put the SMP measurement files (```*.pnt```) into the ```indata``` folder 

3. Run the snowmicrotoolset.py script from the same directory where ```indata``` and ```outdata``` are located. 

From Windows command line:
```
python __init__.py
```

4. All output CSV/PNG files will be sent to the ```outdata``` folder


## License
This dataset is licensed under the [Open Government License of Canada](http://open.canada.ca/en/open-government-licence-canada)
and is subject to the [Copyright Act of Canada](http://laws-lois.justice.gc.ca/eng/acts/C-42/index.html). Additional information can be found at the [Government of Canada's Open Government portal](http://open.canada.ca)
