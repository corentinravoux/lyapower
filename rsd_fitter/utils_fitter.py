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



import logging, time
import numpy as np


#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################

lambdaLy = 1215.673123130217





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




def error_estimator(power,model="uncorrelated",**kwargs):
    if(model == "uncorrelated"):
        return(error_estimator_uncorrelated(power,**kwargs))
    else:
        raise KeyError("model of error estimator not available")


def error_estimator_uncorrelated(power,**kwargs):
    epsilon = return_key(kwargs,"epsilon",None)
    bin_count = return_key(kwargs,"bin_count",0)
    if((bin_count is None)|(epsilon is None)): return KeyError("Need bin_count and epsilon")
    return(power*((1/np.sqrt(bin_count)) + epsilon))
