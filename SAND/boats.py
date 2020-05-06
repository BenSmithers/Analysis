import numpy as np
from math import sqrt
pi = np.pi

"""
Gets an order of magnitude estimate for number of neutrinos on target as a result of other nuclear powered boats in the ocean 
"""

def sci(number):
    order = int(np.log10(number))
    return("{:.2f}e{}".format(number/(10**order),order))

total_reactor_emission = 2.e20 # neuts/second
diff_reactor_emission = total_reactor_emission/(4*pi)

earth_area = 3.61e8 # km^2 
detector_radius = 1e-3 # km
n_boats = 50
unc_boat = 20

#max_distance = sqrt( earth_area/pi )
max_distance = 1300. # km
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
# uncertainty on position will be root(max_distance/12.)
unc = sci(sum(neuts_on_target*((np.sqrt(max_distance/12.)/max_distance) + float(unc_boat)/n_boats ) + sqrt(diff_reactor_emission)/diff_reactor_emission ))

print("{} +/- {}".format(sci(sum(neuts_on_target)), unc))
