try:
    import os
    import struct
    import datetime
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
except ImportError as err:
    print err
    

class SMP(object):
    def __init__(self, binaryfile):
        self.header, self.units = self.retrieve_header(binaryfile)
        self.data = self.extract_data(binaryfile, self.header)
    
    def retrieve_header(self, pnt_file):
        '''
        ripped straight from https://sourceforge.net/projects/pyntreader/
        SMP project site: http://www.slf.ch/ueber/organisation/schnee_permafrost/projekte/SnowMicroPen/index_EN
        '''
        #header construction name, type, start, length, unit
        construct = [
    		['Version','H',0,2, "-"],
    		['Tot Samples','i',2,4, "-"],
    		['Samples Dist [mm]','f',6,4, "mm"],
    		['CNV Force [N/mV]','f',10,4, "N/mV"],
    		['CNV Pressure [N/bar]','f',14,4, "N/bar"],
    		['Offset [N]','H',18,2, "N"],
    		['Year','H',20,2, "y"],
    		['Month','H',22,2, "m"],
    		['Day','H',24,2, "d"],
    		['Hour','H',26,2, "h"],
    		['Min','H',28,2, "min"],
    		['Sec','H',30,2, "s"],
    		['X Coord','d',32,8, "deg"],
    		['Y Coord','d',40,8, "deg"],
    		['Z Coord','d',48,8, "deg"],
    		['Battery [V]','d',56,8, "V"],
    		['Speed [mm/s]','f',64,4, "mm/s"],
    		['Loopsize','l',68,4, "-"],
    		['Waypoints','10l',72,40, "-"],
    		['Calstart','10H',112,20, "-"],
    		['Calend','10H',132,20, "-"],
    		['Length Comment','H',152,2, "-"],
    		['Comment','102s',154,102, "-"],
    		['File Name','8s',256,8, "-"],
    		['Latitude','f',264,4, "deg"],
    		['Longitude','f',268,4, "deg"],
    		['Altitude [cm]','f',272,4, "cm"],
    		['PDOP','f',276,4, "-"],
    		['Northing','c',280,1, "-"],
    		['Easting','c',281,1, "-"],
    		['Num Sats','H',282,2, "-"],
    		['Fix Mode','H',284,2, "-"],
    		['GPS State','c',286,1, "-"],
    		['reserved 1','x',187,1, "-"],
    		['X local','H',288,2, "deg"],
    		['Y local','H',290,2, "deg"],
    		['Z local','H',292,2, "m"],
    		['Theta local','H',294,2, "deg"],
    		['reserved 2','62x',296,62, "-"],
    		['Force Samples','l',358,4, "-"],
    		['Temperature Samples','l',362,4, "-"],
    		['Kistler Range [pC]','H',366,2, "pC"],
    		['Amp Range [pC]','H',368,2, "pC"],
    		['Sensitivity [pC/N]','H',370,2, "pC/N"],
    		['Temp Offset [N]','h',372,2, "Celsius"],
    		['Hand Op','H',374,2, "-"],
    		['Diameter [um]','l',376,4, "um"],
    		['Overload [N]','H',380,2, "N"],
    		['Sensor Type','c',382,1, "-"],
    		['Amp Type','c',383,1, "-"],
    		['SMP Serial','H',384,2, "-"],
    		['Length [mm]','H',386,2, "mm"],
    		['reserved 3','4x',388,4, "-"],
    		['Sensor Serial','20s',392,20, "-"],
    		['Amp Serial','20s',412,20, "-"],
    		['reserved 4 ','80x',432,80, "-"]
        ]
        
        # read in the raw binary data from pnt_file
        with open(pnt_file, 'rb') as in_raw:
            raw = in_raw.read()
        
        # generate list of header values
        values = []
        for f in construct:
            frmt = '>' + f[1]
            start = f[2]
            end = start + f[3]
            try:
                value = struct.unpack(frmt, raw[start:end])[0]
            except:
                value = ""
                pass
            values.append(value)
        
        # generate list of header value names   
        names = [row[0] for row in construct]
        
        # format the names and values as a dict
        header = dict(zip(names, values))
        
        # strip empty-length Comment value strings, or just remove trailing whitespace
        if header['Length Comment'] == 0:
            header['Comment'] = ""
        else:
            header['Comment'] = header['Comment'].split("\x00")[0]
        
        # strip trailing whitespace from some value strings
        header['File Name'] = header['File Name'].split("\x00")[0]	
        header['Amp Serial'] = header['Amp Serial'].split("\x00")[0]
        header['Amp Type'] = header['Amp Type'].split("\x00")[0]
        header['Sensor Serial'] = header['Sensor Serial'].split("\x00")[0]
        header['Sensor Type'] = header['Sensor Type'].split("\x00")[0]
        header['Northing'] = header['Northing'].split("\x00")[0]
        header['Easting'] = header['Easting'].split("\x00")[0]
        header['GPS State'] = header['GPS State'].split("\x00")[0]
        
        # add sign to coordinates if required
        if header["Northing"] == "S":
            header["Latitude"] = -header["Latitude"]
        if header["Easting"] == "W":
            header["Longitude"] = -header["Longitude"]
        
        # measurement units per value
        units = dict(zip(names, [row[4] for row in construct]))
        return header, units
    
    def extract_data(self, pnt_file, header_info):
        '''
        ripped straight from https://sourceforge.net/projects/pyntreader/
        SMP project site: http://www.slf.ch/ueber/organisation/schnee_permafrost/projekte/SnowMicroPen/index_EN
        '''
        # read in the raw binary data from pnt_file
        with open(pnt_file, 'rb') as in_raw:
            raw = in_raw.read()
        try:
            # the starting binary chunk index for the data samples
            start = 512
            # the ending chunk index
            end = header_info['Force Samples'] * 2 + start
            # binary storage format for unpacking
            # big-endian short integers (huh?)
            frmt = '>' + str(header_info['Force Samples']) + 'h'
            data = struct.unpack(frmt, raw[start:end])
        except:
            print "ERROR: unable to read data points from {}".format(os.path.basename(pnt_file))
        else:
            dx = header_info['Samples Dist [mm]']
            data_x = np.arange(0, len(data)) * dx
            data_y = np.asarray(data) * header_info['CNV Force [N/mV]']
            data = np.column_stack([data_x, data_y])
            
            #print "Read {} data points from {}".format(len(data_y), os.path.basename(pnt_file))
            return data
        
    def plot_quicklook(self, out_png):
        '''
        Plot a quick depth/force profile using matplotlib, exporting it to a PNG
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[:,1], self.data[:,0], 'k-', linewidth=0.8)
        x_lim = int(np.ceil(self.data[:,1].max()))
        y_lim = int(np.ceil(self.data[:,0].max()))
        ax.set_ybound(-int(y_lim*.05), y_lim + int(y_lim*.05))
        if x_lim >= 10:
            ax.set_xbound(-int(x_lim*.05), x_lim + int(x_lim*.05))
            ax.set_xticks(range(0, x_lim+1, 5))
            ax.set_xticks(range(0, x_lim+1, 1), minor=True)
        else:
            ax.set_xbound(-int(x_lim*.1), x_lim + int(x_lim*.1))
            ax.set_xticks(range(0, x_lim+1, 1))
        ax.set_xlabel('Force (N)')
        ax.invert_yaxis() # Measures Depth from the top of snowpack, so invert Y axis
        ax.set_ylabel('Depth (mm)')
        ax.set_title(os.path.splitext(os.path.basename(out_png))[0])
        plt.savefig(out_png)
        plt.close(fig)
            
    
    def as_dataframe(self):
        '''
        may be useful at some point in the future
        '''
        return pd.DataFrame(data=self.data, columns=['Depth', 'Force'])
    
    def export_to_csv(self, out_csv):
        '''
        export the depth/force columns to a csv with included 7-line header
        '''
        serial = self.header['SMP Serial']
        sampleTotal = int(self.data.shape[0]) # TOT_SAMPLE field from header is not the number of measurements!
        lat = self.header['Latitude']
        lon = self.header['Longitude']
        tstamp = datetime.datetime(self.header['Year'], self.header['Month'], self.header['Day'], 
                                   self.header['Hour'], self.header['Min'], self.header['Sec'])
        csv_header = " SMP Serial: {}\n {}\n {}\n Lat: {}\n Lon: {}\n Total Samples: {}\n Depth (mm),Force (N)"
        np.savetxt(out_csv, self.data, delimiter=',', comments='#', fmt='%.6f',
                      header=csv_header.format(serial, tstamp.strftime("%Y-%m-%d"), tstamp.strftime("%H:%M:%S"), lat, lon, sampleTotal))
    
    

if __name__ == "__main__":
    workdir = os.path.abspath(os.path.dirname(__file__))
    input_data = os.path.join(workdir, 'indata')
    output_data = os.path.join(workdir, 'outdata')
    pnt_list = []
    # walk through subdirectories of input_data, looking for SMP .pnt files
    for root, folders, files in os.walk(input_data):
        for f in files:
            pnt_list.append(os.path.join(root, f))
    
    for pnt in pnt_list:
        dirname, basename = os.path.split(pnt)
        dirname = dirname.split("\\")
        unique_key = os.path.splitext(basename)[0]
        ################################ WARNING ###############################
        ## The following variables assume a very specific directory structure ##
        ## that was defined by the project lead                               ##
        ################################ WARNING ###############################
        try:
            code = dirname[-3][-1]
            site = dirname[-1]
            date = dirname[-2].split("_")[-1]
        except IndexError:
            print "Error: cannot extract site information from filepath"
            exit(1)
        outcsv = "SMP{}_{}_{}_{}.csv".format(code, site, date, unique_key)
        outcsvabs = os.path.join(output_data, outcsv)
        outpng = outcsvabs.replace(".csv", ".png")
        # still-alive msg
        print outcsv, "{}/{}".format(pnt_list.index(pnt)+1, len(pnt_list))
        p = SMP(pnt)
        if not os.path.isfile(outcsvabs): p.export_to_csv(outcsvabs)
        if not os.path.isfile(outpng): p.plot_quicklook(outpng)
