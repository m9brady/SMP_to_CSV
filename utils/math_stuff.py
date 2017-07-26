import numpy as np
from pylab import rms_flat

def gen_msCoef(dens, ssa):
    '''
    Use the supplied Density and Specific Surface Area params to
    produce bias-minimized microstructure coefficients
    '''
    pass


# Inspired by http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html
def rolling_window(a, window, fun=None, pad=False):
    '''
    Returns moving window bins of size window. Set pad to NaN fill to size of a.
    
    Parameters
    ----------
    a : 
        Input numpy array. Used in SMP class for Force measurement windowing.
    window : {int}
        The size (in bins) of the rolling window.     
    fun :  
        The window function to apply. It currently accepts np.median or np.mean. Default is ``np.mean``.
    pad : {bool}
        Decides whether or not to NaN-pad the outer edges of the rolling window array to match the original input array size. Default is ``False``.
    '''
    window = (np.ceil(window) // 2 * 2 + 1).astype(int) #round up to next odd number
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    rWindow = np.lib.stride_tricks.as_strided(a, shape=shape, strides=a.strides + (a.strides[-1],))
    if fun is not None:
        try:
            rWindow = fun(rWindow, axis=-1)
        except:
            raise Exception('Must provide the proper arguments to the supplied window function: {}'.format(fun.__module__ + " " + fun.__name__))
    if pad:  # slow if fun == None
        padSize = np.abs(rWindow.shape[0] - a.shape[0])/2
        rWindow = np.lib.pad(rWindow, (padSize, padSize), 'constant', constant_values=np.nan)
    return rWindow


# ripped straight from https://sourceforge.net/projects/pyntreader/
# SMP project site: http://www.slf.ch/ueber/organisation/schnee_permafrost/projekte/SnowMicroPen/index_EN
def xcorr(x, y=None, maxlags=None, norm='biased'):
    """Cross-correlation using numpy.correlate
    
    Estimates the cross-correlation (and autocorrelation) sequence of a random
    process of length N. By default, there is no normalisation and the output
    sequence of the cross-correlation has a length 2*N+1. 
    
    :param array x: first data array of length N
    :param array y: second data array of length N. If not specified, computes the 
        autocorrelation. 
    :param int maxlags: compute cross correlation between [-maxlags:maxlags]
        when maxlags is not specified, the range of lags is [-N+1:N-1].
    :param str option: normalisation in ['biased', 'unbiased', None, 'coeff']
     
    The true cross-correlation sequence is
    
    .. math:: r_{xy}[m] = E(x[n+m].y^*[n]) = E(x[n].y^*[n-m])

    However, in practice, only a finite segment of one realization of the 
    infinite-length random process is available.
    
    The correlation is estimated using numpy.correlate(x,y,'full'). 
    Normalisation is handled by this function using the following cases:

        * 'biased': Biased estimate of the cross-correlation function
        * 'unbiased': Unbiased estimate of the cross-correlation function
        * 'coeff': Normalizes the sequence so the autocorrelations at zero 
           lag is 1.0.

    :return:
        * a numpy.array containing the cross-correlation sequence (length 2*N-1)
        * lags vector
        
    .. note:: If x and y are not the same length, the shorter vector is 
        zero-padded to the length of the longer vector.
               
    .. rubric:: Examples
    
    .. doctest::
    
        >>> from spectrum import *
        >>> x = [1,2,3,4,5]
        >>> c, l = xcorr(x,x, maxlags=0, norm='biased')
        >>> c
        array([ 11.])
    
    .. seealso:: :func:`CORRELATION`.  
    """
    N = len(x)
    if y == None:
        y = x
    assert len(x) == len(y), 'x and y must have the same length. Add zeros if needed'
    assert maxlags <= N, 'maxlags must be less than data length'
    
    if maxlags == None:
        maxlags = N-1
        lags = np.arange(0, 2*N-1)
    else:
        assert maxlags < N
        lags = np.arange(N-maxlags-1, N+maxlags)
              
    res = np.correlate(x, y, mode='full')
    
    if norm == 'biased':
        Nf = float(N)
        res = res[lags] / float(N)    # do not use /= !! 
    elif norm == 'unbiased':
        res = res[lags] / (float(N)-abs(np.arange(-N+1, N)))[lags]
    elif norm == 'coeff':        
        Nf = float(N)
        rms = rms_flat(x) * rms_flat(y)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]

    lags = np.arange(-maxlags, maxlags+1)        
    return res, lags
