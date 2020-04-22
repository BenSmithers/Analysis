#!/usr/bin/python3.6
'''
This script plots the fluxes output by the convolve cpp script
'''

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from cross_section_test import get_diff_flux
import os
import nuSQUIDSpy as nsq

# colormap thing
cmap = plt.get_cmap('viridis')
n_colors = 6
def get_color(which, how_many=n_colors):
    return( cmap( float(which)/how_many ) )


# define a couple utility functions
def get_flavor( key ):
    '''
    take a flux dictionary key and return the nusquids flavor type
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

flavors = ['E', 'Mu', 'Tau']
neuts = ['nu', 'nuBar']
currents = ['NC', 'CC']

def get_index( key ):
    '''
    Returns the column in the atmosphere data file for a given key
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))
    split = key.split('_')
    flavor = split[0]
    variety = split[1]

    flav_index = flavors.index(flavor)
    variety_index = neuts.index(variety)
    return( 2 + int( flav_index + len(flavors)*variety_index) )


curr = nsq.NeutrinoCrossSections_Current.NC

# load the data in
print("Extracting Data")
data = np.loadtxt(os.path.join( os.path.dirname(__file__), 'atmosphere.txt'), dtype=float, comments='#',delimiter=' ')
n_energies = 700
n_angles = 100
assert( len(data) == (n_energies*n_angles))

# this funnny indexing is a result of the way I output the data from nuSQuIDS
# it loops through energies for each angle
print("Building Energy and Angle Arrays")
energies = [data[i][0] for i in range(n_energies)]
angles = [data[n_energies*i][1] for i in range(n_angles)]

# let's fill out some flux functions
# in the data file, these data are written in a big list. But that's not a very handy format
# so I'm converting these into 2D arrays
print("Caltulating Fluxes with XS")
fluxes = {}
muon_ones = np.array([ 0. for energy in range(n_energies) ])
not_muon  = np.array([ 0. for energe in range(n_energies) ])
for flav in flavors:
    for neut in neuts:
        for curr in currents:
            key = flav+'_'+neut + '_'+curr
            if flav=='Mu' and curr=='CC':
                continue 

            fluxes[ key ] = [sum([ data[energy+angle*n_energies][get_index(key)] for angle in range(n_angles)])*get_diff_flux(energy, get_flavor(key), get_neut(key), get_curr(key)) for energy in range(n_energies)]
            
            if curr=='NC' and flav=='Mu':
                muon_ones += np.array(fluxes[key])
            else:
                not_muon += np.array(fluxes[key])

scale_e = 10**np.array(energies)

plot_all = False
if plot_all:
    it = 0
    for key in fluxes:
        plt.plot( scale_e, fluxes[key], color=get_color(it),label=key )
        it+=1

    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('reallycool.png', dpi=400)

# plot the stuff!
plt.plot(scale_e, muon_ones, color=get_color(0,2), label="Muon Cascades")
plt.plot(scale_e, not_muon, color=get_color(1,2), label="Other Cascades")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('muon_not_muon.png', dpi=400)
