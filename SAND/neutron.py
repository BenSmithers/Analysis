import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

"""
Works to get an order-of-magnitude approximation for the number of neutrinos-on-target as a result of neutron absorption of Oxygen in water around nuclear reactor
"""

def sci( number ):
    order = int(np.log10(number))
    return("{:.4f}e{}".format( number/(10**order),order))

barn = 10**-28 # meters

# all from PDG 
ncl_water = 83.33/100. #meters
ncl_lead = 10.05/100. # meters
ncl_iron = 10.37/100. # meters 

detector_radius = 1.0 # meters

print("Lead attenuation, 1 meter: {}".format(1.-np.exp(-1./ncl_lead)))
print("Iron attenuation, 1 meter: {}".format(1.-np.exp(-1./ncl_iron)))

baseline = 10. # meters

flux_total = 1.e17 
flux_diff  = flux_total/(8.*4*np.pi*np.pi)

unshielded = 2*np.pi/3.
partial_flux = unshielded*flux_diff 

distances = np.linspace(0,100, 2001)
widths = distances[1:]-distances[:-1]
centers = 0.5*(distances[1:]+distances[:-1])
dist_dist = ncl_water*partial_flux*np.exp( -1*centers/ncl_water )*widths
# now, we need the number that actually interacted with /oxygen/
oxygen_dist = dist_dist * (16./(16+2))

# then this thing emits a neutrino isotropically, so we multiply this by the odds that the emitted neutrino hits 
sr_ratios = (np.pi*(detector_radius**2))/(4*np.pi*np.sqrt( baseline**2 + centers**2))
nus_on_target = oxygen_dist*sr_ratios
print("Estimated nus on target: {}".format(sci(sum(nus_on_target))))


