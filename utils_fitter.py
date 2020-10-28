#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date : 17/05/2019

Author: Corentin Ravoux

Description : Void and Over-density finder for Lya Tomographic maps.
Watershed and Simple Spherical techniques are available.
Tested on irene and cobalt (CCRT)
"""



#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################



import numpy as np
import logging, time
import scipy



#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################

lambdaLy = 1215.673123130217


def bin_ndarray(ndarray, new_shape, operation='mean'):
    """
    From : https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array/29042041
    Bins an ndarray in all axes based on the target shape, by summing or
    averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
    [102 110 118 126 134]
    [182 190 198 206 214]
    [262 270 278 286 294]
    [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg','gauss']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
        elif operation.lower() in ["gauss"]:
            if i!=0 : raise KeyError("gaussian mean is not available for dim higher than 1")
            from scipy import signal
            newndarray = np.zeros(new_shape)
            gaussian_weights = signal.gaussian(int(ndarray.shape[1]),int(ndarray.shape[1])/4)
            for j in range(len(ndarray)):
                newndarray[j] = np.average(ndarray[j],axis=0,weights=gaussian_weights)
            ndarray =newndarray
    return ndarray




def create_log(name="Python_Report",log_level="info"):
    log = Logger(name=name,log_level=log_level)
    log.setup_logging()
    return(log)

def create_report_log(name="Python_Report",log_level="info"):
    log = Logger(name=name,log_level=log_level)
    log.setup_report_logging()
    return(log)

_logging_handler = None

class Logger(object):
    
    def __init__(self,name="Python_Report",log_level="info"):
        self.name = name
        self.log_level = log_level


    def setup_logging(self):
    	"""
        Taken from https://nbodykit.readthedocs.io/
    	Turn on logging, with the specified level.
    	Parameters
    	----------
    	log_level : 'info', 'debug', 'warning'
    		the logging level to set; logging below this level is ignored
    	"""
    
    	# This gives:
    	#
    	# [ 000000.43 ]   0: 06-28 14:49  measurestats	INFO	 Nproc = [2, 1, 1]
    	# [ 000000.43 ]   0: 06-28 14:49  measurestats	INFO	 Rmax = 120
    
    	levels = {"info" : logging.INFO,"debug" : logging.DEBUG,"warning" : logging.WARNING}
    
    	logger = logging.getLogger();
    	t0 = time.time()
    
    
    	class Formatter(logging.Formatter):
    		def format(self, record):
    			s1 = ('[ %09.2f ]: ' % (time.time() - t0))
    			return s1 + logging.Formatter.format(self, record)
    
    	fmt = Formatter(fmt='%(asctime)s %(name)-15s %(levelname)-8s %(message)s',
    					datefmt='%m-%d %H:%M ')
    
    	global _logging_handler
    	if _logging_handler is None:
    		_logging_handler = logging.StreamHandler()
    		logger.addHandler(_logging_handler)
    
    	_logging_handler.setFormatter(fmt)
    	logger.setLevel(levels[self.log_level])



    def setup_report_logging(self):
        levels = {"info" : logging.INFO,"debug" : logging.DEBUG,"warning" : logging.WARNING}
        logging.basicConfig(filename=self.name, filemode='w',level=levels[self.log_level],format='%(asctime)s :: %(levelname)s :: %(message)s')



    @staticmethod
    def add(line,level="info"):
        if(level=="info"):
            logging.info(line)
        if(level=="warning"):
            logging.warning(line)
        if(level=="debug"):
            logging.debug(line)

    @staticmethod
    def close():
        logging.shutdown()






def latex_float(float_input,decimals_input="{0:.2g}"):
    """
    example use: 
    import matplotlib.pyplot as plt
    plt.figure(),plt.clf()
    plt.plot(np.array([1,2.]),'ko-',label="$P_0="+latex_float(7.63e-5)+'$'),
    plt.legend() 
    """
    float_str = decimals_input.format(float_input)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


def return_key(dictionary,string,default_value):
    return(dictionary[string] if string in dictionary.keys() else default_value)



