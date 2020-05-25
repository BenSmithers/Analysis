import nuSQUIDSpy as nsq
from math import sqrt

import numpy as np

"""
This defines a few utility functions for my plotting script.

I moved this over here so the main plotter wasn't too busy
"""

def get_width( which_list ):
    """
    Takes a list 'which_list' of floats of length N, considered the centers of some bins
    Returns a length N list of floats for the widths of the bins these centers correspond to

    Arg 'which_list' must be monotonically increasing or decreasing
    """
    if not (isinstance(which_list, list) or isinstance(which_list, np.ndarray) or isinstance(which_list, tuple)):
        raise TypeError("Expected list-like, got {}".format(type(which_list)))

    # cast as a list
    if not isinstance(which_list, list):
        use = list(which_list)
    else:
        use = which_list

    n_bins = len(which_list)
    if n_bins<=1:
        raise ValueError("'which_list' should be longer than length {}".format(len(which_list)))

    increasing = which_list[1] > which_list[0]
    for i in range(n_bins-1):
        if not(increasing == (which_list[i+1] > which_list[i])):
            raise TypeError("Arg 'which_list' should be monotonically increasing or decreasing")

    n_bins = len(which_list)
    widths = [ 0. for i in range(n_bins)]
    for i in range(n_bins):
        if i==0:
            widths[i] = abs(use[1]-use[0])
        elif i==(n_bins-1):
            widths[i] = abs(use[-1]-use[-2])
        else:
            widths[i] = abs(0.5*(use[i+1]-use[i-1]))
    return(np.array(widths))

def get_exp_std( widths, probabilities, values ):
    """
    This takes a discretely binned probability density function representing a true continuous one,
    integrates it to get the expectation value and standard deviation of the distribution

    If the probability density isn't normalized, this function will normalize it

    TODO: get asymmetric error bars! 

    RETURNS; means, y+ error, y- error 
    """
    if not (len(widths)==len(values) or len(values)==len(prob_density)):
        raise ValueError("The args should have the same lengths")

    norm = sum(probabilities*widths)
    prob_density = [ value/norm for value in probabilities]

    mean=0.
    for item in range(len(widths)):
        # P(x)*X
        mean+= (widths[item]*prob_density[item])*values[item]

    var=0.
    for item in range(len(widths)):
        var+= (widths[item]*prob_density[item])*(mean-values[item])**2
    sigma = sqrt(var)

    return( mean, sigma, sigma)

# define a couple utility functions
def get_flavor( key ):
    '''
    take a flux dictionary key and return the nusquids flavor type
    
    The dictionary key will be like "electorn_stuff_stuff"
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))

    part = key.split('_')[0].lower()
    if part in ['e', 'eleectron']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.electron )
    elif part in ['mu', 'muon']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.muon )
    elif part in ['tau']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.tau )
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))

def get_neut( key ):
    '''
    Takes a flux dictionary key and returns the nusquids neutrino type 
    (anti neutrino or vanilla neutrino)
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))   
    part = key.split('_')[1].lower()

    if part in ['nu', 'matter','neutrino']:
        return(nsq.NeutrinoCrossSections_NeutrinoType.neutrino)
    elif part in ['nubar', 'antimatter', 'antineutrino']:
        return(nsq.NeutrinoCrossSections_NeutrinoType.antineutrino)
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))

def get_curr(key):
    '''
    Takes a flux dictionary key and returns the nusquids neutrino type 
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))   
    part = key.split('_')[2].lower()

    if part in ['neutral', 'nc']:
        return(nsq.NeutrinoCrossSections_Current.NC)
    elif part in ['charged', 'cc']:
        return(nsq.NeutrinoCrossSections_Current.CC)
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))


