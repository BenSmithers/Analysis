from math import exp, sqrt, pi
import numpy as np
import os

"""
This script is here to approximate the uncertainties in going from "energy deposited" to "energy reconstructed"
"""

rtwo = 1./sqrt(2*pi)

datafile = "reconstruction.txt"
data = np.loadtxt(os.path.join(os.path.dirname(__file__), datafile), dtype=float, delimiter=",")
data = data.transpose()
#data[1] = data[1]*data[0] 
# data[0] is energy
# data[1] is sigma

from utils import get_closest


def get_odds(deposited, reconstructed):
    """
    Takes an energy deposited and energy reconstructed.
    Loads the datafile and reads off the uncertainty from the second column. 
    The data is in %E_depo, so we scale this to be a ratio and then by that deposted energy
    """
    if not isinstance(deposited, (float,int)):
        raise Exception()
    if not isinstance(reconstructed, (float,int)):
        raise Exception()

    sigma = get_closest( energy, data[0], data[1])*deposited*0.01
    # now, we assume that the uncertainty follows a Gaussian distribution, and calculate the PDF here
    prob = rtwo*(1./sigma)*exp(-0.5*((deposited-reconstructed)/sigma)**2)
    
    return(prob)


