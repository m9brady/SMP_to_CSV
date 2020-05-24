import numpy as np
from scipy import signal


'''
Helper functions to work with the SMP data

'''

def xcorr_microstructure(smp, obs):
    '''
    Most of the SMP profiles we collect are close, but not directly through where the SSA samples are colelcted
    There is therefore some uncertainty in the equivalent heights between two sensors.
    This function uses a cros-correlation to extimate an offset between
    
    @param smp: An SMP class object
    @param obs: Co-located SSA estimates as 2D numpy array where obs[:,0] are the estimate height
        refrenced to ground/ice and obs[:,1] are the SSA estimates/observations
    return: Offset in mm to maximize correlation between the two datasets
    
    '''
    if smp.shotNoise is None:
            return "SMP object with shotNoise estimates required"
     
    obsHeight = obs[:,0]
    obsSSA = obs[:,1]
    
    #Need to flip the SSA depth values to snow height
    #If the ground position is wrong, this is wrong.
    smpHeight = (max(smp.shotNoise[:,5])-smp.shotNoise[:,5])/10
    obsSSAEqv = np.empty(smpHeight.size)
    obsSSAEqv[:] = np.NAN

    for i_step in range(0,obsSSA.size):
        idxMatch = min(range(len(smpHeight)), key=lambda i: abs(smpHeight[i]-obsHeight[i_step]))
        obsSSAEqv[idxMatch] = obsSSA[i_step]
        
    #Find the NaN values and interplate between the SSA obs
    nans, x= np.isnan(obsSSAEqv), lambda z: z.nonzero()[0]
    obsSSAEqv[nans]= np.interp(x(nans), x(~nans), obsSSAEqv[~nans])  
    
    obsDetrend = signal.detrend(obsSSAEqv - obsSSAEqv.mean())
    smpDetrend = signal.detrend(smp.microStructure[:,3] - smp.microStructure[:,3].mean())
    
    xCorrSSA = np.correlate(smpDetrend, obsDetrend, mode='full')
    xCorrMaxIdx = np.argmax(xCorrSSA)-obsDetrend.size
    return(xCorrMaxIdx*1.25)
