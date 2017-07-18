import numpy as np
# Inspired by http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
def rolling_window(a, window, fun=None, pad=False):
    '''
    Returns moving window bins of size window. Set pad to NaN fill to size of a.
    
    Parameters
    ----------
    a : 
        input numpy array.
    window : {int}
        the size (in bins) of the rolling window.     
    fun :  
        the window function to apply. It currently accepts np.median or np.mean. Default is ``np.mean``.
    pad : {bool}
        decides whether or not to NaN-pad the outer edges of the rolling window array to match the original input array size. Default is ``False``.
    '''
    window = (np.ceil(window) // 2 * 2 + 1).astype(int) #round up to next odd number
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    rWindow = np.lib.stride_tricks.as_strided(a, shape=shape, strides=a.strides + (a.strides[-1],))
    if fun is not None:
        try:
            rWindow = fun(rWindow, -1)
        except:
            raise Exception('Must provide the proper arguments to the supplied window function: {}'.format(fun.__module__ + " " + fun.__name__))
    if pad:  # slow if fun == None
        padSize = np.abs(rWindow.shape[0] - a.shape[0])/2
        rWindow = np.lib.pad(rWindow, (padSize, padSize), 'constant', constant_values=np.nan)
    return rWindow


def gen_msCoef(dens, ssa):
    '''
    Use the supplied Density and Specific Surface Area params to
    produce bias-minimized microstructure coefficients
    '''
    pass
