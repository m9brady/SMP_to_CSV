# SMP_to_CSV
Converts header and data files ("\*_Header.txt" and "\*_Data.txt") into CSV containing all data points (Depth and Pen Force) plus a multi-line header containing time, latitude and longitude information.

## Requirements
*Assumes you have Python 2.7 installed with pandas and matplotlib modules*

1. Create two folders "indata" and "outdata" in the same directory as *SMP_to_CSV.py* 

2. Put the header ("\*_Header.txt") and data ("\*_Data.dat") files into the "indata" folder

3. Run the SMP_to_CSV.py script from the same directory where "indata" and "outdata" are located. 

From Windows command line:
```
python SMP_to_CSV.py
```

4. All output CSV files will be sent to the "outdata" folder
