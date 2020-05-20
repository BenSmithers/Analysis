#!/usr/bin/python3.6
'''
This script plots the fluxes output by the convolve cpp script
'''

import numpy as np
import matplotlib
# Need to use agg since Tk isn't on the cobalts??? 
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
    variety = split[1] # nu or nu-bar 

    flav_index = flavors.index(flavor)
    variety_index = neuts.index(variety)
    return( 2 + int( flav_index + len(flavors)*variety_index) )



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
an_width = [0. for i in range(n_angles)]
for i in range(n_angles):
    if i ==0:
        an_width[0] = abs( np.arccos(angles[1]) - np.arccos(angles[0]))
    else:
        an_width[i] = abs( np.arccos(angles[i]) - np.arccos(angles[i-1]))

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
                # skip tracks 
                continue 
            scale_with = 1.
#            scale_with = get_diff_flux((10**energies[energy]), get_flavor(key), get_neut(key), get_curr(key))

            fluxes[ key ] = [sum([ data[energy+angle*n_energies][get_index(key)]*an_width[angle] for angle in range(n_angles)])*scale_with*(en_width[energy])*2*np.pi for energy in range(n_energies)]
            
            if curr=='NC' and flav=='Mu':
                muon_ones += np.array(fluxes[key])
            else:
                not_muon += np.array(fluxes[key])

def get_flux( energy, key, use_overflow = False):
    '''
    interpolates between entries in the flux dictionary to return the flux at arbitrary energy
    interpolation done in log-space of energy! 


    returns DOUBLE  (0.0 if beyond scope of data)
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
    CC_E_energy = 0.98 # ratio of energy recovered in a CC electron cascade
    event = 100*const.GeV # observed energy of cascade 
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
    print("In Energy Increasing" if in_energies[1]>in_energies[0] else "In Energy Decreasing")
    from_muon = np.array([ 0. for energy in range(n_bins) ])
    from_not = np.array([ 0. for energy in range(n_bins) ])
    from_diffy = {}
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                    continue 
                   
                # in-energy needs to be changed. 
                # if it's a normal hadronic cascade (muon NC)

                # if it's an electron cascade, then that means this total energy is split between the hadronic and electric component... uugh. Then we need to consider the energy of the electron
                # OH! Then we just care about the total cross section! 

                key = flav+'_'+neut + '_'+curr
                if flav=='E' and curr=='CC':
                    # in this case, we recover the energy from the entire cascade, so we essentially exactly know the energy of the parent neutrino
                    # and so get the total cross-section/flux from that energy and put it all in one bin
                    force_initial_energy = event/CC_E_energy 
                    from_diffy[key] = [0.0 for j in range(n_bins)]
                    which_bin = n_bins - 1
                    while in_energies[which_bin]<force_initial_energy:
                        which_bin-=1
                    from_diffy[key][which_bin] = get_flux(np.log10(in_energies[which_bin]), key)*get_diff_flux(in_energies[which_bin], get_flavor(key), get_neut(key), get_curr(key))*widths[which_bin] 
                else:
                    from_diffy[key] = [get_flux(np.log10(in_energies[j]), key)*get_diff_flux(in_energies[j], get_flavor(key), get_neut(key), get_curr(key),in_energies[j]-event,0.)*widths[j] for j in range(n_bins)]
                # this evaluates in (N / GeV)
                if curr=='NC' and flav=='Mu':
                    from_muon += np.array(from_diffy[key])
                else:
                    from_not += np.array(from_diffy[key])

    
    isolation = False
    if isolation:
        norm = sum(from_muon)+sum(from_not)
        # these are properly normalized so that you can integrate over part of the trend to get the probability the event was of that type 
        for i in range(n_bins):
            from_muon[i]=from_muon[i]/(widths[i]*norm)
            from_not[i] =from_not[i]/(widths[i]*norm)


        plt.plot(in_energies/const.GeV, from_muon, color=get_color(0,2),label="Muon Origin")
        plt.plot(in_energies/const.GeV, from_not, color=get_color(1,2),label="Other Origin")
    else:
        norm = sum([ sum(from_diffy[key]) for key in from_diffy ])
        n_colors = len(list(from_diffy.keys()))
        counter = 0
        for key in from_diffy:
            from_diffy[key] = [value/norm for value in from_diffy[key]]
            plt.plot( in_energies/const.GeV, from_diffy[key], color=get_color(counter, n_colors),label=key)
            counter+=1

    plt.xscale('log')
    plt.yscale('log')
    plt.title("{:.2f}GeV Cascade Rates".format(event/const.GeV))
    plt.xlabel("Parent Neutrino Energy [GeV]")
    plt.ylabel(r"Probability [GeV$^{-1}$]")
    plt.legend()
    plt.savefig("wow.png",dpi=400)
