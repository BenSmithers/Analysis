import numpy as np
from math import sqrt
pi = np.pi

"""
Gets an order of magnitude estimate for number of neutrinos on target as a result of other nuclear powered boats in the ocean 
"""

def sci(number):
    order = int(np.log10(number))
    return("{:.4f}e{}".format(number/(10**order),order))

total_reactor_emission = 2.e20 # neuts/second
diff_reactor_emission = total_reactor_emission/(4*pi)

earth_area = 3.61e8 # km^2 
detector_radius = 1e-3 # km
n_boats = 200

#max_distance = sqrt( earth_area/pi )
max_distance = 3500. # km
print("Effective Max Radius: {}".format(max_distance))

def prob(r):
    """
    returns in m^-1

    this should be multiplied by bin width! 
    """
    return( 1.5*(max_distance**(-3./2))*np.sqrt(r) )


dist_bins = np.linspace(0, max_distance, 2001)
widths = dist_bins[1:]-dist_bins[:-1]
centers=0.5*(dist_bins[1:]+dist_bins[:-1])

probs = prob(centers)
neuts_on_target = n_boats*total_reactor_emission*probs*(pi*detector_radius**2)/(4*pi*centers**2)
print("{}".format(sci(sum(neuts_on_target))))
