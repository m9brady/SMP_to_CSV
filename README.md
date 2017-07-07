# SMP_to_CSV
Converts SMP measurement files ("*.pnt") into CSV containing all data points (Depth and Penetration Force) plus a multi-line header containing time, latitude, longitude and other site-specific information. A depth/force profile quicklook is also generated as a PNG file.

## Disclaimer
The SMP class methods "retrieve_header" and "extract_data" are from the excellent [SnowMicroPyn](https://sourceforge.net/projects/pyntreader/files/) tool by Sascha Grimm. SnowMicroPyn is licensed under the [GNU General Public License version 3.0 (GPLv3)](https://sourceforge.net/directory/license:gplv3/)

The "detect_peaks" and associated "_plot" helper functions are from the [Biomechanics and Motor Control](https://github.com/demotu/BMC/blob/master/functions/detect_peaks.py) repo by Marcos Duarte. BMC is licensed under the [MIT Licence](https://github.com/demotu/BMC/blob/master/LICENSE.txt)

## Requirements and Usage
*Assumes you have Python 2.7 installed with numpy, scipy, pandas and matplotlib modules*

### Option A: Running the tool standalone

1. Create two folders ```indata``` and ```outdata``` in the same directory as SMP_to_CSV.py

2. Put the SMP measurement files ("*.pnt") into the ```indata``` folder 

3. Run the SMP_to_CSV.py script from the same directory where ```indata``` and ```outdata``` are located. 

From Windows command line:
```
python SMP_to_CSV.py
```

4. All output CSV/PNG files will be sent to the ```outdata``` folder

### Option B: Importing the SMP class definition

In your Python interactive session or script tool, first ensure that the path to ```SMP_to_CSV.py``` is properly appended to your system path. To use the SMP class in your work, use the following command:

```
from SMP_to_CSV import SMP
```

## License
This dataset is licensed under the [Open Government License of Canada](http://open.canada.ca/en/open-government-licence-canada)
and is subject to the [Copyright Act of Canada](http://laws-lois.justice.gc.ca/eng/acts/C-42/index.html). Additional information can be found at the [Government of Canada's Open Government portal](http://open.canada.ca)
