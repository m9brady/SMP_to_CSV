import os
import datetime
from glob import glob

import pandas as pd
import matplotlib.pyplot as plt

# generate a depth/force profile for the current SMP site
def save_quicklook(data_vals, out_png):
    d_list = []
    for d in data_vals:
        d_list.append([float(x.strip()) for x in d])
    dframe = pd.DataFrame(data=d_list, columns=['Depth', 'Force'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    dframe.plot(x='Force', y='Depth',linewidth=0.8, color='k', legend=False, ax=ax)
    x_lim = int(pd.np.ceil(dframe['Force'].max()))
    ax.set_xbound(-1, x_lim)
    ax.set_xticks(range(0, x_lim, 5))
    ax.set_xticks(range(0, x_lim, 1), minor=True)
    ax.set_xlabel('Force (N)')
    ax.set_ylabel('Depth (mm)')
    ax.set_title(os.path.splitext(os.path.basename(out_png))[0])
    fig.savefig(out_png)
    

# use the current shell directory
workdir = os.getcwd()

input_data = os.path.join(workdir, 'indata')
output_data = os.path.join(workdir, 'outdata')

# assumes that files are alphabetically sorted
in_hdr_files = sorted(glob(os.path.join(input_data, '*_Header.txt')))
in_dat_files = sorted(glob(os.path.join(input_data, '*_Data.dat')))

if not len(in_hdr_files) == len(in_dat_files):
    print "ERROR: Mismatched number of data and header files!\n.dat count: {}\nheader count: {}".format(len(in_dat_files), len(in_hdr_files))
    exit(1)

# pairs header and dat files together
#pair_list = zip(in_hdr_files, in_dat_files) #old method, doesn't ensure matches between header and data files
pair_list = []
for hdr_file in in_hdr_files:
    pair = [hdr_file, hdr_file.replace("_Header.txt", "_Data.dat")]
    if not os.path.isfile(pair[1]):
        print "Cannot locate data file associated with {}".format(os.path.basename(hdr_file))
    else:
        pair_list.append(pair)

for file_pair in pair_list:
    # check for existing csv file
    csv_file = os.path.join(output_data, os.path.basename(file_pair[1]).replace('_Data.dat', '.csv'))
    if os.path.isfile(csv_file):
        print "output file exists for data file: {}".format(os.path.basename(file_pair[1]))
        continue
    
    # first, read the header for relevant info
    with open(file_pair[0], 'rb') as in_hdr:
        for row in in_hdr:
            # Time Information (no timezone info, assume UTC?)
            if row.startswith('Day'):
                day = int(row.split()[-1])
            elif row.startswith('Month'):
                month = int(row.split()[-1])
            elif row.startswith('Year'):
                year = int(row.split()[-1])
            elif row.startswith('Hour'):
                hour = int(row.split()[-1])
            elif row.startswith('Min'):
                minute = int(row.split()[-1])
            elif row.startswith('Sec'):
                second = int(row.split()[-1])
            # LatLon info
            elif row.startswith('Latitude'):
                lat = float(row.split()[-1])
            elif row.startswith('Longitude'):
                lon = float(row.split()[-1])
        # generate datetime object from time info
        tstamp = datetime.datetime(year, month, day, hour, minute, second)
    
    # second, read the datafile for the observations
    data_list = [] # for holding our SMP data rows
    with open(file_pair[1], 'rb') as in_dat:
        # skip header row
        in_dat.next()
        for row in in_dat:
            depth, force = row.split()
            data_list.append([depth, force+'\n'])
            
    # third, write to csv
    with open(csv_file, 'wb') as out_csv:
        # write the multi-line header first
        out_csv.write("{}: {} \nLatitude: {}\nLongitude: {}\n{},{}\n".format('Timestamp', tstamp.isoformat(), lat, lon, 'Depth (mm)', 'Force (N)'))
        # write the SMP data rows
        for row in data_list:
            out_csv.write(",".join(row))
            
    print "{} generated".format(os.path.basename(csv_file))
    save_quicklook(data_list, csv_file.replace(".csv",".png"))