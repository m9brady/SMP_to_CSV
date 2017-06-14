try:
    import os
    import struct
    import datetime
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from tools_ext import detect_peaks
    from scipy import signal 
except ImportError as err:
    print err
    

class SMP(object):
    def __init__(self, binaryfile):
        self.header, self.units = self.retrieve_header(binaryfile)
        self.data = self.extract_data(binaryfile, self.header)
        self.subset = self.data[self.pick_surf():,:]
        #self.shotNoise = self.shot_noise()
           
    def pick_surf(self, sWindow = 1000):
        smooth_data = self.moving_average(self.data[:,1],sWindow)
        noiseMed = np.median(smooth_data[sWindow:sWindow*2]) # using median for now
        ind = detect_peaks(smooth_data, mph=noiseMed*2,mpd=sWindow, show=True)
        return(ind[0])
    
    def moving_average(self, a, n=1000):
        '''
        Smooths the SMP force data with a running average with window size of n
        Pads the front of the array with nan values to return  an array of the same length
        '''
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ret[n - 1:] = ret[n - 1:] / n
        ret[:n - 1] = np.nan
        return ret
    
    def est_microstucture(self, coef):
        '''
        Estimate the microstucture properties based on Proksch et al 2015
        Ported from various sources written by Martin Proksch and Josh King
        '''
        if self.shotNoise is None:
            return "Shot noise parameters are required"
        densSmp = np.empty(self.shotNoise[:,0].size, dtype=np.float64)
        longCorSmpEx = np.empty(self.shotNoise[:,0].size, dtype=np.float64)
        longCorSmpC = np.empty(self.shotNoise[:,0].size, dtype=np.float64)
        ssaSmp = np.empty(self.shotNoise[:,0].size, dtype=np.float64)
        for i in range(0,self.shotNoise[:,0].size-1):
           densSmp[i] = coef['a1'] + (coef['a2']*np.log(self.shotNoise[i,0])) + (coef['a3']*np.log(self.shotNoise[i,0])*self.shotNoise[i,4]) + (msCoef['a4']*self.shotNoise[i,4])
           phiSmp = densSmp[i]/ 916.7
           longCorSmpEx[i] = coef['b1'] + (coef['b2'] * self.shotNoise[i,4]) + (coef['b3'] * np.log(self.shotNoise[i,0]))
           longCorSmpC[i] = coef['c1'] + (coef['c2'] * self.shotNoise[i,4]) + (coef['c3'] * np.log(self.shotNoise[i,0]))
           ssaSmp[i] = (4*(1-phiSmp)) / longCorSmpC[i]
        return(np.stack((densSmp,longCorSmpEx,longCorSmpC,ssaSmp),axis=-1))
    
    def shot_noise(self, A_cone = 19.6, window_size_mm = 2, overlap_mm = 0.5):
        '''
        Estimate the shot noise paramters based on Löwe and van Herwijnen, 2012
        Ported from various sources written by Martin Proksch, J-B Madore, and Josh King
        '''
        windowSize = int(round(window_size_mm/self.header['Samples Dist [mm]']))
        #stepSizeMM = overlap_mm*window_size_mm
        stepSize = int(overlap_mm*windowSize)
        #nWindows=int(np.floor(self.subset[:,1].size/windowSize))
        nSteps=int(np.floor(self.subset[:,1].size/stepSize-1))
        
        z = np.empty(nSteps, dtype=np.float64)
        medf_z = np.empty(nSteps, dtype=np.float64)
        delta = np.empty(nSteps, dtype=np.float64)
        lam = np.empty(nSteps, dtype=np.float64)
        f_0 = np.empty(nSteps, dtype=np.float64)
        L = np.empty(nSteps, dtype=np.float64)
        #plt.plot(p.data[:,0],p.data[:,1])

        for i_step in range(1,nSteps):
            z_min=(i_step-1)*stepSize+1
            z_max=(i_step-1)*stepSize+windowSize
            f_z=self.subset[z_min:z_max,1]
            N=len(f_z);
        
            # calc z-vector:
            z[i_step] = (i_step-1)*stepSize + stepSize
            z[i_step] = round(z[i_step]*p.header['Samples Dist [mm]']*100)/100
          
            # calc median penetration force:
            medf_z[i_step] = np.median(f_z)
            
            # calc shot noise
            c1 = np.mean(f_z)
            c2 = np.var(f_z) * (len(f_z) - 1) / len(f_z)
            A = signal.detrend(f_z-c1)
            C_f = np.correlate(A, A, mode='full')
            
            #Normalize xcorr by n-lags for 'unbiased'
            maxlag=N-1
            lags = np.append(np.linspace(N-maxlag,N-1,N-1),N)
            lags = np.append(lags,np.linspace(N-1,N-maxlag,N-1))
            lags = lags *np.repeat(1, C_f.size)
            C_f = C_f/lags #normalize by n-lag
            
            #Shot noise parameters
            delta[i_step] = -3/2 * C_f[N] / (C_f[N+1] - C_f[N]) * self.header['Samples Dist [mm]']; # eq. 11 in Löwe and van Herwijnen, 2012  
            lam[i_step] = 4/3 * np.power(c1,2) / c2 / delta[i_step] # eq. 12 in Löwe and van Herwijnen, 2012
            f_0[i_step] = 3/2 * c2 / c1  # eq. 12 in Löwe and van Herwijnen, 2012
            L[i_step] = np.power(A_cone/lam[i_step],1./3.) 
            
            self.shotNoise = np.stack((medf_z,delta,lam,f_0,L),axis=-1)
            
        return(self.shotNoise)
    
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
    os.chdir(workdir)
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
        
        #placeholder microstucture coefficients
        msCoef = {'a1': 420.47, #+-8.31 kg/m-3
            'a2': 102.47, #+-4.24 N-1
            'a3': -121.15, #+-10.65 N-1mm-1
            'a4': -169.96, #+-18.70 mm-1
            'b1': 0.0715, #+-0.0058 mm
            'b2': 0.299, #+-0.011 mm-1
            'b3': 0.0149, #+-0.0018 N-1
            'c1': 0.131, #+-0.0081 mm
            'c2': 0.155, #+-0.015 mm-1
            'c3': 0.0291} #+-0.0024 N-1
        #shot_noise(self)
        #est_microstucture(self, msCoef)