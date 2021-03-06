# -*- coding: utf-8 -*-
from __future__ import print_function  # py2/3 print compatibility

import datetime
import os
import struct

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal, stats

from .utils import detect_peaks, rolling_window


class SMP(object):
    '''
    Placeholder for description of SMP class
    '''
    def __init__(self, binaryfile):
        self.header, self.units = self.retrieve_header(binaryfile)
        self.data = self.extract_data(binaryfile, self.header)
        self.qFlags = {
            'snowSurfaceFound': False, # by default, we have not defined the surfaces within the SMP measurement data
            'soilSurfaceFound': False, 
            'shortRun': False, # too-shallow SMP measurement
            'airShot': False, # air-shot refers to SMP device test without any snow measurements
            # SMP signal quality classes as per [Pielmeier & Marshall 2009] http://www.sciencedirect.com/science/article/pii/S0165232X0900113X
            'C1': False, # C1 flag, No detected errors
            'C2': False, # C2 flag, Trend or offset in absolute SMP force
            'C3': False, # C3 flag, Dampened or disturbed SMP force micro-variance
            'C4': False
        } # C4 flag, both C2 and C3
        self.shotNoise = None # may not always need to estimate shot noise, so we initialize as None to save processing time
        self.microStructure = None # may not always need to estimate microstructure, so we initialize as None to save processing time
        self.subset = None
        self.snowSurf = self.data[0, 0] # by default, first depth value
        self.soilSurf = self.data[-1, 0] # by default, last depth value
    
    
    def __str__(self):
        '''
        string representation of SMP object : 
            site name [time of measurement] (number of raw measurements)
            
        Example
        -------
        >>> x = SMP('S34M1358.pnt')
        >>> print(x)
        >>> 'S34M1358 [2017-02-20T19:30:33] (290400 raw measurements)'
        '''
        d = datetime.datetime(
            self.header['Year'], self.header['Month'], 
            self.header['Day'], self.header['Hour'], 
            self.header['Min'], self.header['Sec']
        )
        return "%s [%s] (%d raw measurements)" % (
            self.header['File Name'], 
            d.isoformat(), 
            self.data.shape[0]
        )


    def as_dataframe(self, use_raw=False):
        '''
        may be useful at some point in the future.
        
        Parameters
        ----------
        use_raw : {bool}    (default=False)
            use the raw SMP data rather than the valid-data subset.
        '''
        if use_raw:
            return pd.DataFrame(data=self.data, columns=['Depth', 'Force'])
        else:
            if self.subset is None:
                print("Data not preprocessed, use filter_raw() method")
                return None
            else:
                return pd.DataFrame(data=self.subset, columns=['Depth', 'Force'])
        
    
    def est_microstructure(self, coef):
        u'''
        Estimate the microstructure properties based on [1]_.
        Ported from various sources written by Martin Proksch and Josh King
        
        References
        ----------
        .. [1] Proksch, M., Löwe, H. and Schneebeli, M. (2015), `Density, specific surface area, and correlation length of snow measured by high-resolution penetrometry`_. J. Geophys. Res. Earth Surf., 120: 346362.
        
        .. _Density, specific surface area, and correlation length of snow measured by high-resolution penetrometry:
            http://dx.doi.org/10.1002/2014JF003266
		'''
        if self.shotNoise is None:
            return "Shot noise parameters are required."
        arrLen = self.shotNoise.shape[0]
        #TODO: Possible refactor to a single 2D array with columns for each property
        densSmp = np.empty(arrLen)      # Density (kg/m^3)
        longCorSmpEx = np.empty(arrLen) # Exponential Correlation Length (mm)
        longCorSmpC = np.empty(arrLen)  # Correlation Length (mm)
        ssaSmp = np.empty(arrLen)       # Specific Surface Area (mm^-1)
        
        log_medf_z = np.log(self.shotNoise[:,0]) # log of the window median force
        L = self.shotNoise[:,4]
        densSmp = coef['a1'] + (coef['a2'] * log_medf_z) + (coef['a3'] * log_medf_z * L) + (coef['a4'] * L)
        phiSmp = densSmp / 916.7
        longCorSmpEx = coef['b1'] + (coef['b2'] * L) + (coef['b3'] * log_medf_z)
        longCorSmpC = coef['c1'] + (coef['c2'] * L) + (coef['c3'] * log_medf_z)
        ssaSmp = (4 * (1 - phiSmp)) / longCorSmpC
            
        self.microStructure = np.column_stack((densSmp, longCorSmpEx, longCorSmpC, ssaSmp))
        return 0 # success
    
    
    def est_shot_noise(self, A_cone=19.6, window_size_mm=2.0, overlap=0.5):
        u'''
        Estimate the shot noise parameters based on [1]_.
        Ported from various sources written by Martin Proksch, J-B Madore, and Josh King.
        
        Parameters
        ----------
        A_cone : {float}    (default=19.6)
            The projected cone area in square-millimeters.
        window_size_mm : {float}    (default=2.0)
            The size of the rolling window in millimeters.
        overlap : {float}   (default=0.5) 
            Used with ``window_size_mm`` to define the step size.
            
        Returns
        -------
        Nx6 numpy array (dtype = float64) where N is the number of windows based on the number of raw observations and the ``window_size_mm``/``overlap`` args.
        
        References
        ----------
        .. [1] Löwe, H., and van Herwijnen, A. (2012), `A Poisson shot noise model for micro-penetration of snow`_, Cold Regions Science and Technology, Volume 70, 2012, Pages 62-70, ISSN 0165-232X. 
        
        .. _A Poisson shot noise model for micro-penetration of snow:
            http://www.sciencedirect.com/science/article/pii/S0165232X11001832 
        '''
        if self.subset is None:
            print("Data not preprocessed, use filter_raw() method")
            return 1 # nonzero return value implies bad stuff went down
        samplesDist = self.header['Samples Dist [mm]']
        windowSize = int(round(window_size_mm / samplesDist))
        #stepSizeMM = overlap_mm*window_size_mm
        stepSize = int(overlap * windowSize)
        #nWindows=int(np.floor(self.subset[:,1].size/windowSize))
        nSteps = int(np.floor(self.subset.shape[0] / stepSize - 1))
        
        #TODO: Possible refactor to a single 2D array with columns for each parameter
        z = np.empty(nSteps)        # Depth-vector (?)
        medf_z = np.empty(nSteps)   # Median Penetration Force
        delta = np.empty(nSteps)    # Intensity of Poisson process(?)
        lam = np.empty(nSteps)      # Rupture strength(?)
        f_0 = np.empty(nSteps)      # Deflection at rupture(?)
        L = np.empty(nSteps)        # Element size(?)
        #plt.plot(p.data[:,0],p.data[:,1])

        for i_step in xrange(nSteps):
            z_min = i_step * stepSize
            z_max = i_step * stepSize + windowSize
            f_z = self.subset[z_min:z_max,1]
            N = len(f_z)
        
            # calc z-vector 
            # TODO: Check this, not sure if its valid
            z[i_step] = round((i_step * stepSize + stepSize) * samplesDist * 100) / 100
          
            # calc median penetration force
            medf_z[i_step] = np.median(f_z)
            
            # calc shot noise
            c1 = f_z.mean()
            c2 = f_z.var()  # var is population var by default with np
            A = signal.detrend(f_z - c1)
            C_f = np.correlate(A, A, mode='full')
            
            # Normalize xcorr by n-lags for 'unbiased'
            maxlag = N - 1
            lags = np.append(np.linspace(N - maxlag, N - 1, N - 1), N)
            lags = np.append(lags, np.linspace(N - 1, N - maxlag, N - 1))
            lags *= np.repeat(1, C_f.size) # MB: lags is just being multiplied by 1? candidate for removal?
            C_f /= lags # normalize by n-lag
            
            # Shot noise parameters
            delta[i_step] = -3. / 2 * C_f[N-1] / (C_f[N] - C_f[N-1]) * samplesDist # eq. 11 in Löwe and van Herwijnen, 2012  
            lam[i_step] = 4. / 3 * c1 ** 2 / c2 / delta[i_step] # eq. 12 in Löwe and van Herwijnen, 2012
            f_0[i_step] = 3. / 2 * c2 / c1  # eq. 12 in Löwe and van Herwijnen, 2012
            L[i_step] = (A_cone / lam[i_step]) ** (1. / 3) 
        
        self.shotNoise = np.column_stack((medf_z, delta, lam, f_0, L, z))
        return 0 # success
    
    
    def export_to_csv(self, outCsv):
        '''
        using the raw data array, export the depth/force columns to a csv with included 7-line header
        '''
        serial = self.header['SMP Serial']
        sampleTotal = int(self.data.shape[0]) # "Tot Samples" field from header is not the number of measurements!
        lat = round(self.header['Latitude'], 10)
        lon = round(self.header['Longitude'], 10)
        time_src = 'GPS' if self.header['GPS State'] == '1' else 'UTC'
        tstamp = datetime.datetime(
            self.header['Year'], self.header['Month'], 
            self.header['Day'], self.header['Hour'], 
            self.header['Min'], self.header['Sec']
        )
        csv_header = " SMP Serial: {}\n Date: {}\n Time ({}): {}\n Lat: {}\n Lon: {}\n Total Samples: {}\n Depth (mm),Force (N)"
        np.savetxt(
            outCsv, self.data, delimiter=',', comments='#', fmt='%.6f',
            header=csv_header.format(serial, tstamp.strftime("%Y-%m-%d"), 
            time_src, tstamp.strftime("%H:%M:%S"), lat, lon, sampleTotal)
        )
    
    
    def extract_data(self, pnt_file, header_info):
        '''
        ripped straight from https://sourceforge.net/projects/pyntreader/
        SMP project site: https://github.com/slf-dot-ch/snowmicropyn
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
            print("ERROR: unable to read data points from {}".format(os.path.basename(pnt_file)))
        else:
            dx = header_info['Samples Dist [mm]']
            data_x = np.arange(0, len(data)) * dx
            data_y = np.asarray(data) * header_info['CNV Force [N/mV]']
            data = np.column_stack([data_x, data_y])
            #debug
            #print "Read {} data points from {}".format(len(data_y), os.path.basename(pnt_file))
            return data
    
    
    def filter_raw(self, zCor=-1):
        '''
        Some of the SMP units are producing negative force values
        This causes havoc when estimating microstructure params
        Give the user the option to replace them with 0s or interpolate from the nearest non-zero values
        Default decision is to prompt the user for case-by-case input (zCor=-1)
        '''
        filteredArr = self.data.copy()
        
        # pass 1: short run (less than 100mm?)
        if filteredArr[-1, 0] < 100:
            print('***Warning: Short run encountered in dataset: {src}'.format(src=self.header['File Name']))
            self.qFlags['shortRun'] = True
            
        # pass 2: air shot
        force_variance = filteredArr[:,1].var(dtype=np.float64, ddof=1) # ddof=1 "...provides an unbiased estimator of the variance of a hypothetical infinite population"
        if force_variance < 0.01:
            print('***Warning: dataset is likely an air-shot: {src}'.format(src=self.header['File Name']))
            self.qFlags['airShot'] = True
                
        # pass 3: negative values
        if any(filteredArr[:,1] < 0):
            if zCor == -1:
                try:
                    zCor = int(raw_input("Negative force values detected. Replace with 0s (1), interpolate (2): "))
                except (SyntaxError, ValueError): # if we can't cast the input string to an int just go with the default
                    print("Invalid entry. Negative forces left uncorrected")
                    return filteredArr
            if zCor == 1:
                print("Negative forces replaced with zeroes for dataset: {src}".format(src=self.header['File Name']))
                filteredArr[filteredArr[:,1] < 0,1] = 0
            elif zCor == 2:
                print("Negative forces being interpolated with valid-data for dataset: {src}".format(src=self.header['File Name']))
                #Adapted from https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
                filteredArr[filteredArr[:,1] < 0,1] = np.nan
                nans, x = np.isnan(filteredArr[:,1]), lambda z: z.nonzero()[0]
                filteredArr[nans,1] = np.interp(x(nans), x(~nans), filteredArr[~nans,1])
            else:
                print("Negative forces left uncorrected") 
        return filteredArr
    
    
    def mask_ice_lenses(self):
        '''
        Minimize the effect of ice lense presence in the profile. Ice lenses are identified by thin sections with
        unusually-sharp increases in the required penetration force.
        
        Operates upon the filtered subset of the raw SMP data.
        '''
        pass
    
    
    #TODO: TESTING!; Use full data record or subset? If subset, check for it
    def outliers(self, windowMM=5, threshold=1, pad=False):
        '''
        Robust Z-Score outliers; can be used to detect ice features and or soil
        
        Parameters
        ----------
        windowMM : {float}    (default=5.0)
            The window size in millimeters.
        threshold : {int}    (default=1)
            For determining whether or not an outlier is present in the dataset.
        pad : {bool}    (default=False)
            Decides whether or not to NaN-pad the resultant array to match the input array (currently self.subset).
        '''
        sWindow = windowMM / self.header['Samples Dist [mm]']
        sWindow = (np.ceil(sWindow) // 2 * 2 + 1).astype(int)
        sBins = rolling_window(self.subset[:,1], sWindow, None, pad=False)
        medFilter = rolling_window(self.subset[:,1], sWindow, fun=np.median, pad=False)
        diff = np.sqrt(np.sum((sBins - medFilter.reshape(medFilter.size, 1))**2, axis=-1))
        mAD = np.nanmedian(diff) # Median Absolute Deviation
        if pad:
            padSize = np.absolute(diff.size - self.subset[:,1].size) / 2
            diff = np.lib.pad(diff, (padSize,padSize), 'constant', constant_values=np.nan)
        rZScore = 0.6745 * diff / mAD # Robust Z-Score
        #This will throw warnings if the pad is applied, equality does not work for NaN
        return rZScore > np.nanmedian(rZScore) + (threshold * np.nanstd(rZScore))
    
    
    def pick_surf(self, surfType, sWindow=1000):
        '''
        To avoid spurious SMP force measurements, this method identifies the 
        beginning(end) of valid SMP force measurements at the snow(soil) 
        surface and returns the depth value in millimeters.
        
        Parameters
        ----------
        surfType : {'snow', 'soil'} 
            Either 'snow' or 'soil' for the respective surface of interest.
        sWindow : {int}     (default=1000)
            The size of the smoothing window for the noise-reducing moving-average approach in bins.
            
        Notes
        -----
        This method should only ever be called *once* per surface type for a given SMP file, since it queries and then subsequently is used to alter the ``self.subset`` SMP attribute.
        '''
        # snow surface
        if surfType == 'snow':
            # first, smooth the filtered data to reduce noise
            smoothedForce = rolling_window(self.subset[:,1], sWindow, fun=np.mean, pad=True)
            noiseMed = np.median(smoothedForce[sWindow:sWindow * 2]) # using median for now
            # identify the pen-force peaks, toggle show to True to display inline plot
            ind = detect_peaks(smoothedForce, mph=noiseMed * 2, mpd=sWindow, show=False)
            # return the **first** pen-force peak (i.e. the top of the snowpack)
            try:
                surfIdx = ind[0]
            # if we can't find any peaks, just return None which gets handled later
            except IndexError:
                print("***Error: Cannot locate the snow surface")
                return self.snowSurf # default value has been setup in __init__ as the first depth reading
            # if all is well, trip the quality flag for determining the snow surface
            self.qFlags['snowSurfaceFound'] = True
            snowSurf = self.subset[surfIdx, 0] # yank the depth value using our new index
            return snowSurf
        # soil surface
        elif surfType == 'soil':
            # adapted from https://sourceforge.net/projects/pyntreader/
            depth = self.subset[:,0]
            force = self.subset[:,1]
            overload = self.header['Overload [N]']
            surfIdx = -1
            soilSurf = depth[surfIdx] # by default, have the soil surface be the deepest available value
            # if the overload value is tripped anywhere in the force data,
            # the SMP may have hit the ground so we do some checks
            if force.max() >= overload:
                surfIdx = np.argmax(force) # index of the maximum force value
                i_thresh = np.where(depth >= depth[surfIdx] - 20)[0][0]
                f_mean = force[:i_thresh].mean()
                f_std = force[:i_thresh].std()
                thresh = f_mean + 5 * f_std
                while force[surfIdx] > thresh:
                    surfIdx -= 10
                # if the current depth is shallower than the max depth,
                # we trip the quality flag for identifying the soil surface
                if depth[surfIdx] < soilSurf:
                    self.qFlags['soilSurfaceFound'] = True
                soilSurf = depth[surfIdx]
            # return depth value only
            return soilSurf
        else:
            raise ValueError("Invalid surfType argument. Must be one of ['snow', 'soil']")   
            
            
    def plot_ms(self):
        '''
        debug-usage, plot the microstructure parameters
        '''
        fig = plt.figure(figsize=(11,8))
        fig.suptitle(self.header['File Name'], fontsize=16)
        ax0 = fig.add_subplot(221, title='Density (kg m$^{-3}$)')
        ax1 = fig.add_subplot(222, title='Exponential Correlation Length (mm)')
        ax2 = fig.add_subplot(223, title='Correlation Length (mm)')
        ax3 = fig.add_subplot(224, title='Specific Surface Area (mm$^{-1}$)')
        ax0.plot(self.microStructure[:,0])
        ax1.plot(self.microStructure[:,1])
        ax2.plot(self.microStructure[:,2])
        ax3.plot(self.microStructure[:,3])
        #fig.tight_layout()
        plt.show()
    
    
    def plot_quicklook(self, outPng):
        '''
        Create a quick depth/force profile of the raw SMP data using matplotlib, exporting it to a PNG
        
        Parameters
        ----------
        outPng : 
            the absolute path to the output PNG file. Uses the current working directory if abspath is not passed.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data[:,1], self.data[:,0], 'k-', linewidth=0.8)
        x_lim = int(np.ceil(self.data[:,1].max()))
        y_lim = int(np.ceil(self.data[:,0].max()))
        ax.set_ybound(-int(y_lim * .05), y_lim + int(y_lim * .05))
        if x_lim >= 10:
            ax.set_xbound(-int(x_lim * .05), x_lim + int(x_lim * .05))
            ax.set_xticks(range(0, x_lim + 1, 5))
            ax.set_xticks(range(0, x_lim + 1, 1), minor=True)
        else:
            ax.set_xbound(-int(x_lim * .1), x_lim + int(x_lim * .1))
            ax.set_xticks(range(0, x_lim + 1, 1))
        ax.set_xlabel('Force (N)')
        ax.invert_yaxis() # Measures Depth from the top of snowpack, so invert Y axis
        ax.set_ylabel('Depth (mm)')
        ax.set_title(os.path.splitext(os.path.basename(outPng))[0])
        if not outPng.endswith('.png'): outPng += ".png" # force it!
        plt.savefig(outPng)
        plt.close(fig)


    def plot_self(self, filterSize=5, plotSurfs=True):
        '''
        debug-usage, plot the raw data, overlaid with subset, overlaid with ``rolling_mean`` function
        with vertical dashed lines denoting detected snow/soil surfaces
        
        **filterSize**: for controlling the mean filter windowsize (Default: 5mm)
        
        **plotSurfs**: decide if detected surfaces should be plotted, if available (Default: True)
        '''
        if self.subset is None:
            sub = self.data[:,1].copy()
        else:
            sub = self.subset[:,1].copy()
            sub = np.concatenate([np.empty(self.data.shape[0] - self.subset.shape[0]), sub])
        sub[np.where(self.data[:,0] < self.snowSurf)[0]] = np.NaN
        e = rolling_window(sub, int(filterSize/self.header['Samples Dist [mm]']), np.mean, pad=True)
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)

        ax.set_ylim(0, self.data[:,0].max())
        ax.plot(self.data[:,1], self.data[:,0], color='red', label='Raw Data')
        ax.plot(sub, self.data[:,0], color='cyan', label='Filtered Subset of Raw', linewidth=2.2)
        ax.plot(e, self.data[:,0], color='darkgreen', label='~{}mm window mean'.format(filterSize))
        if plotSurfs:
            if self.qFlags['snowSurfaceFound']:
                ax.axhline(self.snowSurf, label='Snow Surface (~{}mm)'.format(round(self.snowSurf,2)), linestyle='dashed', color='k', linewidth=0.8)
                ax.text(ax.get_xlim()[1]*.2, self.snowSurf-5, '\nSnow Surface')
            if self.qFlags['soilSurfaceFound']:
                ax.axhline(self.soilSurf, label='Estimated Soil Surface (~{}mm)'.format(round(self.soilSurf,2)), linestyle='dashed', color='k', linewidth=0.8)
                ax.text(ax.get_xlim()[1]*.2, self.soilSurf-5, '\nSoil Surface')
            
        ax.set_ylabel('Depth (mm)')
        ax.set_xlabel('Force (N)')
        ax.invert_yaxis()
        ax.legend(framealpha=1)
        ax.set_title(self.header['File Name'], fontsize=16)
        fig.tight_layout()
        plt.show()
    

    #  **WORK IN PROGRESS - PROBABLY WRONG**
    def quality_pass(self):
        '''
        Using the raw SMP data, assign values to self.qFlags signal quality classes as per [1]_.
        
        References
        ----------
        .. [1] Pielmeier, C., and Marshall, H-P. (2009), `Rutschblock-scale snowpack stability derived from multiple quality-controlled SnowMicroPen measurements`_, Cold Regions Science and Technology, Volume 59, Issue 2, 2009, Pages 178-184, ISSN 0165-232X
        
        .. _Rutschblock-scale snowpack stability derived from multiple quality-controlled SnowMicroPen measurements:
            http://www.sciencedirect.com/science/article/pii/S0165232X0900113X
        
        '''
        f = self.data[:,1].copy()
        d = self.data[:,0].copy()
        _,_,r,_,_ = stats.linregress(f,d)
        
        # first, check for presence of linear trend (indicates a systematic error with the SMP device)
        if ((r**2) >= 0.7) or (self.header['Offset [N]'] != 0):  # also check for offset recorded by SMP device
            self.qFlags['C2'] = True
            
        # second, check for "dampened or disturbed SMP force micro-variance"
        window = int(np.floor(5 / self.header['Samples Dist [mm]']))
        maxDeviation = 0.02 # threshold
        rollingVariance = rolling_window(f, window, np.var, pad=False)
        if np.any(rollingVariance > maxDeviation):
            self.qFlags['C3'] = True
            
        # C4 flag is triggered if both C2 and C3 are True
        self.qFlags['C4'] = True if (self.qFlags['C2'] and self.qFlags['C3']) else False
        
        # C1 flag is True if the other flags are false
        self.qFlags['C1'] = True if not(self.qFlags['C2'] or self.qFlags['C3']) else False
        return 0 # success
    
    
    def retrieve_header(self, pnt_file):
        '''
        ripped straight from https://sourceforge.net/projects/pyntreader/
        SMP project site: https://github.com/slf-dot-ch/snowmicropyn
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
        
        # try to convert bytes to str where possible so python3 is ok
        for k,v in header.items():
            if type(v) == bytes:
                try:
                    header[k] = v.decode('utf-8')
                except ValueError:
                    pass

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
