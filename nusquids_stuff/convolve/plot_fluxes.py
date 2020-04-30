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
import sys

const = nsq.Const()

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
if energies[1]>energies[0]:
    print("growing")
else:
    print("decreasing")
en_width = [0. for i in range(n_energies)]
for i in range(n_energies):
    if i==0:
        en_width[i] = (10**energies[1]- 10**energies[0])/const.GeV
    else:
        en_width[i] = (10**energies[i] - 10**energies[i-1])/const.GeV
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
            scale_with = 1.
#            scale_with = get_diff_flux((10**energies[energy]), get_flavor(key), get_neut(key), get_curr(key))

            fluxes[ key ] = [sum([ data[energy+angle*n_energies][get_index(key)] for angle in range(n_angles)])*scale_with*(en_width[energy])*4*np.pi for energy in range(n_energies)]
            
            if curr=='NC' and flav=='Mu':
                muon_ones += np.array(fluxes[key])
            else:
                not_muon += np.array(fluxes[key])

def get_flux( energy, key, use_overflow = False):
    '''
    interpolates between entries in the flux dictionary to return the flux at arbitrary energy
    interpolation done in log-space of energy! 


    returns NONE if energy is beyond edges of flux data
    returns DOUBLE otherwise 
    '''
    if not (key in fluxes):
        raise ValueError("Bad key {}".format(key))
    if not (isinstance(energy, float) or isinstance(energy, int)):
        raise TypeError("Expected {}, not {}".format(float, type(energy)))

    

    # check if it's outside the extents
    if energy < energies[0]:
        if use_overflow:
            return(energies[0])
        else:
            return(0)
    if energy > energies[-1]:
        if use_overflow:
            return(energies[n_energies - 1])
        else:
            return(0)
    
    # should execute in O(N) time 
    upper_boundary = 1
    while energy>energies[upper_boundary]:
        upper_boundary += 1
    lower_boundary = upper_boundary - 1

    # sanity check... 
    # essentially makes sure that the energies are monotonically increasing 
    if not ((energies[lower_boundary] <= energy) and (energies[upper_boundary] >= energy)):
        print("energy: {}".format(energy))
        print("lower bound: {}".format(energies[lower_boundary]))
        print("upper bound: {}".format(energies[upper_boundary]))
        print("indices: {}, {}".format(lower_boundary, upper_boundary))
        raise Exception()

    # logarithmic interpolation
    # linear in log-energy space 
    y2 = fluxes[key][upper_boundary]
    y1 = fluxes[key][lower_boundary]
    x2 = energies[upper_boundary]
    x1 = energies[lower_boundary]
    
    flux_value = energy*((y2-y1)/(x2-x1) ) + y2 -x2*((y2-y1)/(x2-x1))
    return(flux_value)

scale_e = 10**np.array(energies)

# this block here is supposed to just plot all the raw fluxes 
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

plt.clf()
# this block here was supposed to compare muon cascade fluxes with non-muon cascade fluxes
casc_compare = False
if casc_compare:
    plt.plot(scale_e, muon_ones, color=get_color(0,2), label="Muon Cascades")
    plt.plot(scale_e, not_muon, color=get_color(1,2), label="Other Cascades")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Incoming Neut Energy [GeV]")
    plt.title("Event Rate for atm flux 10GeV->1PeV")
    plt.ylabel("Event Rate [s$^{-1}$]")
    plt.legend()
    plt.savefig('muon_not_muon.png', dpi=400)
plt.clf()


# now I'm trying tofix what the final state energy is, and look at the relative rates of initial neutrino energy that would give rise to this event 
cool_plot = True
if cool_plot:
    event = 100*const.GeV
    by_min = -5
    by_max = 0
    n_bins = 200
    # the initial neutrino energy we are considering
    in_energies = np.array([event/(1.-y) for y in (1-np.logspace(by_min, by_max, n_bins))])
    widths = [0. for i in range(n_bins)]
    for i in range(n_bins):
        if i==0:
            widths[i]=abs(in_energies[1]-in_energies[0])/const.GeV
        else:
            widths[i]=abs(in_energies[i]-in_energies[i-1])/const.GeV

    from_muon = np.array([ 0. for energy in range(n_bins) ])
    from_not = np.array([ 0. for energy in range(n_bins) ])
    from_diffy = {}
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if flav=='Mu' and curr=='CC':
                    continue 
                   

                key = flav+'_'+neut + '_'+curr
                from_diffy[key] = [get_flux(np.log10(in_energies[j]), key)*get_diff_flux(in_energies[j], get_flavor(key), get_neut(key), get_curr(key), event,0.)*widths[j] for j in range(n_bins)]
                # this evaluates in (N / GeV)
                if curr=='NC' and flav=='Mu':
                    from_muon += np.array(from_diffy[key])
                else:
                    from_not += np.array(from_diffy[key])
    norm = sum(from_muon)+sum(from_not)
    # these are properly normalized so that you can integrate over part of the trend to get the probability the event was of that type 
    for i in range(n_bins):
        from_muon[i]=from_muon[i]/(widths[i]*norm)
        from_not[i] =from_not[i]/(widths[i]*norm)


    plt.plot(in_energies/const.GeV, from_muon, color=get_color(0,2),label="Muon Origin")
    plt.plot(in_energies/const.GeV, from_not, color=get_color(1,2),label="Other Origin")
    plt.xscale('log')
    plt.yscale('log')
    plt.title("{:.2f}GeV Cascade Rates".format(event/const.GeV))
    plt.xlabel("Parent Neutrino Energy [GeV]")
    plt.ylabel(r"Probability [GeV$^{-1}$]")
    plt.legend()
    plt.savefig("wow.png",dpi=400)
