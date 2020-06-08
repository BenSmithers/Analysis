#!/usr/bin/python3.6
'''
This script plots the fluxes output by the convolve cpp script

But how does it work?
 1. Raw flux data from the included mceq+nuSQuIDS flux is loaded in by the Data Object. This object has functionality for sampling from the flux at arbitrary parent neutrino energy by interpolating between nsq points.
 2. This flux data is convolved with differential cross sections to create a incoming V outgoing doubly differential flux-rate array
 3. The DD flux array is converted to be a incoming vs cascade doubly differential flux-rate array, respecting cascade energy differences wrt flavor/interaction

The heavy lifting goes on in the generate_doubly_diff_fluxes function 

Different "modes" are used to make different plots. 
'''

from optparse import OptionParser
import sys

parser = OptionParser()
mode_str = "0 - plot muon vs not muon\n1 - plot all the keys\n 2 - plot 2D hist of parent vs cascade"
parser.add_option("-m", "--mode",
                dest="mode",
                default=0,
                type=int,
                help=mode_str)
parser.add_option("-l", "--load_stored",
                dest="load_stored",
                default=False,
                action="store_true",
                help="Should I try to load stored data rather than regenerate it?")
parser.add_option("-d", "--debug", 
                dest="debug",
                default=False,
                action="store_true",
                help="activates debug options")
options, args = parser.parse_args()
mode = options.mode
load_stored = options.load_stored
debug = options.debug

do_norm=False # deprecated 

recognized_modes = [0,1,2,3,4,5,6,7,-1,8,9]
if mode not in recognized_modes: 
    raise ValueError("Unrecognized Mode: {}".format(mode))
if mode in [-1,0,1,2,4,7]:
    raise DeprecationWarning("Mode {} is deprecated".format(mode))


'''
Modes
    0 - plot muon vs not muon (normed)
    1 - plot all the keys (normed)
    2 - plot 2D hist of parent vs cascade 
    4 - PLOT 2d HIST OF parent vs cascade, doubly differential
    5 - The plot showing most probable energy for muon events and non-muon events 
    6 - The plot just showing must probable energy for all of them
'''

print("Configuration...")
print("    In Mode {}".format(mode))
print("    Will Normalize" if do_norm else "    No Norm - true fluxes")
print("    Will Load Data" if load_stored else "    Will Generate Data")

# data analysis
import numpy as np

# file system, control
import os
import pickle
from warnings import warn

#plotting imports
import matplotlib
# Need to use agg since Tk isn't on the cobalts??? 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker #used for log-scale contourplots 

from cross_section_test import get_diff_xs
import nuSQUIDSpy as nsq

# specialty-made utility functions
from utils import get_flavor, get_neut, get_curr, get_exp_std, get_width

const = nsq.Const()
# colormap thing
cmap = plt.get_cmap('coolwarm')
n_colors = 6
def get_color(which, how_many=n_colors):
    return( cmap( float(which)/how_many ) )

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
    flavor = split[0] # E Mu Tau
    variety = split[1] # nu or nu-bar 

    flav_index = flavors.index(flavor) # 0, 1, or 2
    variety_index = neuts.index(variety) # 0 or 1
    return( 2 + int( flav_index + len(flavors)*variety_index) )

class Data:
    """
    This is used as a container for the data loaded in from nuSQuIDS. 
    The main benefit of this is that the objects used by the interpolator are kept in a sterile scope, 
        and we don't have to worry about accidentally renaming an important object!  
    
    It loads it up into a convenient format for access, and provides a function for interpolating what is loaded. 
    """
    def __init__(self, filename='atmosphere.txt'):
        """
        Loads in the specified nuSQuIDS datafile. 

        Creates a "flux" dictionary for each type of neutrino and interaction. This is in units of N/s/GeV/cm2
        """
        print("Extracting Data")
        data = np.loadtxt(os.path.join( os.path.dirname(__file__), 'atmosphere.txt'), dtype=float, comments='#',delimiter=' ')
        n_energies = 700
        n_angles = 100
        assert( len(data) == (n_energies*n_angles))

        # this funny indexing is a result of the way I output the data from nuSQuIDS
        # it loops through energies for each angle
        print("Building Flux and Energy Arrays")
        self.energies = [10**data[i][0] for i in range(n_energies)]
        self.growing = self.energies[1]>self.energies[0]
        print("Growing" if self.growing else "Decreasing")
        en_width = get_width(self.energies)/const.GeV
        angles = [data[n_energies*i][1] for i in range(n_angles)]
        an_width = get_width(np.arccos(angles))

        # let's fill out some flux functions
        # in the data file, these data are written in a big list. But that's not a very handy format
        # so I'm converting these into 2D arrays
        self.fluxes = {}
        for flav in flavors:
            for neut in neuts:
                for curr in currents:
                    key = flav+'_'+neut + '_'+curr
                    if flav=='Mu' and curr=='CC':
                        # skip tracks 
                        continue 
                    
                    #self.fluxes[ key ] = [sum([ data[energy+angle*n_energies][get_index(key)]*an_width[angle] for angle in range(n_angles)])*(en_width[energy])*2*np.pi for energy in range(n_energies)]
                    self.fluxes[ key ] = [sum([ data[energy+angle*n_energies][get_index(key)]*an_width[angle] for angle in range(n_angles)])*2*np.pi for energy in range(n_energies)]

    def get_flux(self, energy, key, use_overflow = False):
        '''
        interpolates between entries in the flux dictionary to return the flux at arbitrary energy
        Energy should be in units of eV
        Flux is in units of /cm2/GeV/s 

        returns DOUBLE  (0.0 if beyond scope of data)
        '''
        if not (key in self.fluxes):
            raise ValueError("Bad key {}".format(key))
        if not (isinstance(energy, float) or isinstance(energy, int)):
            raise TypeError("Expected {}, not {}".format(float, type(energy)))

        

        # check if it's outside the extents
        if energy < self.energies[0]:
            if use_overflow:
                return(self.energies[0])
            else:
                return(0)
        if energy > self.energies[-1]:
            if use_overflow:
                return(self.energies[n_energies - 1])
            else:
                return(0)
        
        # should execute in O(N) time 
        upper_boundary = 1
        while energy>(self.energies[upper_boundary]):
            upper_boundary += 1
        lower_boundary = upper_boundary - 1

        # sanity check... 
        # essentially makes sure that the energies are monotonically increasing 
        if not ((self.energies[lower_boundary] <= energy) and (self.energies[upper_boundary] >= energy)):
            print("energy: {}".format(energy))
            print("lower bound: {}".format(self.energies[lower_boundary]))
            print("upper bound: {}".format(self.energies[upper_boundary]))
            print("indices: {}, {}".format(lower_boundary, upper_boundary))
            raise Exception()
        
        # linear interpolation 
        y2 = self.fluxes[key][upper_boundary]
        y1 = self.fluxes[key][lower_boundary]
        x2 = self.energies[upper_boundary]
        x1 = self.energies[lower_boundary]
        slope = (y2-y1)/(x2-x1)

        flux_value = energy*slope + y2 -x2*slope
        return(flux_value)

data = Data()

scale_e = np.array(data.energies)

# this block here is supposed to just plot all the raw fluxes 
if debug:
    print("Making Debug plot")
    binny = len(scale_e)
    taus = np.zeros(binny)
    muon = np.zeros(binny)
    ele = np.zeros(binny)
    
    for key in data.fluxes:
        flav=str(key).split('_')[0]
        if 'Tau'==flav:
            taus+=data.fluxes[key]
        elif "E"==flav:
            ele+=data.fluxes[key]
        elif "Mu"==flav:
            muon+=data.fluxes[key]
        else:
            raise Exception("Ya Done goofed")


    plt.plot( scale_e, muon, color=get_color(0,3),label="Muons")
    plt.plot(scale_e, taus, color=get_color(1,3),label="Taus")
    plt.plot(scale_e, ele, color=get_color(2,3), label="Electrons")
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('raw fluxes.png', dpi=400)
    print("saved raw fluxes")
    plt.close()

n_bins = 400


def get_distribs_for_cascade_energy(cascade_energy, event_widths, event_energies):
    """
    Returns an un-normalized dictionary of fluxes for each of the progenitor particles.
    The thing returns in units of /GeV/s (energy of incoming particle)

    Param "cascade_energy" refers to the measured energy of the cascade in eV. FLOAT 
    Param "event_width" refers to the width of the bin for the cascade_energy's energy. Needed to properly scale rates. [GEV]
    Param "event_energies" is the list of possible progenitor neutrino energies to consider. in eV. LIST
    """
    raise DeprecationWarning("This function is deprecated")

    if not (isinstance(cascade_energy, float) or isinstance(cascade_energy,int)):
        raise TypeError("Expected number-like, not {}".format(type(cascade_energy)))
    if not (isinstance(event_energies, list) or isinstance(event_energies, np.ndarray) or isinstance(event_energies, tuple)):
        raise TypeError("Expected list-like, not {}".format(type(event_energies)))
    rising = event_energies[1] > event_energies[0]
    print("Looking at {:.2f} GeV cascade".format(cascade_energy/const.GeV))
    print("Event energies range from {} GeV to {} GeV".format(min(event_energies)/const.GeV,max(event_energies)/const.GeV))
    CC_E_energy = 0.98 # ratio of energy recovered in a CC electron cascade
    widths = get_width(event_energies)/const.GeV
    from_diffy = {}
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                    continue                    
                # in-energy needs to be changed. 
                # if it's a normal hadronic cascade (muon NC)
                key = flav+'_'+neut + '_'+curr
                from_diffy[key] = np.zeros(n_bins)
                if curr=='CC':
                    if flav=='E':
                        print("Injecting Spike")
                        #inject single spike
                        forced_initial_energy = cascade_energy / CC_E_energy
                        which_bin = 0 if rising else n_bins-1
                        while event_energies[which_bin]<forced_initial_energy:
                            which_bin = which_bin+1 if rising else which_bin-1
                            if (which_bin==n_bins or which_bin==-1):
                                break                
                        if (which_bin!=n_bins and which_bin>=0):
                            # in this case, there is only one parent particle energy which could cause this cascade. So we get the total cross section!
                            # N/(cm2 s GeV )     cm2
                            from_diffy[key][which_bin] = data.get_flux(event_energies[which_bin], key)*get_diff_xs(event_energies[which_bin], get_flavor(key), get_neut(key), get_curr(key))
                    else: # this is a track
                        continue
                else: # cool! 
                    lepton_energy = (event_energies - cascade_energy)
                    print("Lepton Energy ranges from {:.2f} to {:.2f}. It is {}".format(min(lepton_energy)/const.GeV, max(lepton_energy)/const.GeV, "increasing" if lepton_energy[0]<lepton_energy[1] else "decreasing"))
                    for j in range(n_bins): # iterate over possible event energies 
                        if lepton_energy[j]>0.:
                            from_diffy[key][j] =data.get_flux(event_energies[j], key)*get_diff_xs(event_energies[j], get_flavor(key), get_neut(key), get_curr(key),lepton_energy[j],0.5)*event_widths[j]
    return(from_diffy, widths)

if mode==-1:
    """
    This mode was made to plot the separate fluxes from the get_distribs_for_cascade_energy function

    That function is deprectated now, and so this one is too. 
    """
    plt.figure(2)
    e_min=50*const.GeV
    e_max=10*const.PeV
    event_energies = np.logspace(np.log10(e_min),np.log10(e_max), n_bins)

    cascade_energy = 100*const.GeV
    lepton_energies = event_energies - cascade_energy
    lepton_widths = get_width(lepton_energies)

    from_diffy, widths = get_distribs_for_cascade_energy(cascade_energy, lepton_widths/const.GeV, event_energies)
    from_elec = np.array([0. for ebin in range(n_bins)])
    from_muon = np.array([0. for ebin in range(n_bins)])
    from_tau  = np.array([0. for ebin in range(n_bins)])
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if (flav=='Mu' or flav=='Tau') and curr=='CC':
                    continue
                key = flav+'_'+neut+'_'+curr
                if flav=='E':
                    from_elec += from_diffy[key]
                elif flav=='Mu':
                    from_muon+=from_diffy[key]
                elif flav=='Tau':
                    from_tau+=from_diffy[key]
                else:
                    raise Exception("Impossible Logic! Flavor {}".format(flav))

    plt.clf()
    plt.plot(event_energies/const.GeV, from_elec, color=get_color(0,3), label="Electrons")
    plt.plot(event_energies/const.GeV, from_muon, color=get_color(1,3), label="Muons")
    plt.plot(event_energies/const.GeV, from_tau , color=get_color(2,3), label="Taus")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Event Energy [GeV]")
    plt.ylabel("Flux by Flavor")
    plt.title("Fluxes leading to a {:.2f} GeV Cascade".format(cascade_energy/const.GeV))
    filename = "Flavor_Debug.png"
    plt.savefig(filename, dpi=400)
    print("Saved {}".format(filename))

if mode==0 or mode==1:
    """
    This was used to make plots showing the probability density of event energy using a fixed cascade energy.

    This is deprecated. 
    """

    e_min = 50*const.GeV
    e_max = 10*const.PeV
    in_energies = np.logspace(np.log10(e_min), np.log10(e_max), n_bins)
    # the initial neutrino energy we are considering

    energy = 100*const.GeV
    casc_widths = get_width(in_energies - energy)
    from_diffy, widths  = get_distribs_for_cascade_energy(energy, casc_widths/const.GeV, in_energies)
    from_muon = np.array([ 0. for ienergy in range(n_bins) ])
    from_not = np.array([ 0. for ienergy in range(n_bins) ])
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                    continue 
                key = flav+'_'+neut + '_'+curr
                if curr=="NC" and flav=="Mu":
                    from_muon += from_diffy[key]
                else:
                    from_not += from_diffy[key]
    if mode==0:
        norm = sum(widths*from_muon)+sum(widths*from_not)
        # these are properly normalized so that you can integrate over part of the trend to get the probability the event was of that type 
        from_muon=(from_muon)/norm
        from_not =(from_not)/norm

        plt.clf()
        plt.plot(in_energies/const.GeV, from_muon, color=get_color(0,2),label="Muon Origin")
        plt.plot(in_energies/const.GeV, from_not, color=get_color(1,2),label="Other Origin")
        
        print("Total Muon probability: {}".format( sum( from_muon*widths )))
        print("Total Not Muon probability: {}".format(sum(from_not*widths)))

    elif mode==1:
        norm = sum([ sum(from_diffy[key]) for key in from_diffy ])
        n_colors = len(list(from_diffy.keys()))
        counter = 0
        plt.clf()
        for key in from_diffy:
            from_diffy[key] = [value/norm for value in from_diffy[key]]
            plt.plot( in_energies/const.GeV, from_diffy[key], color=get_color(counter, n_colors),label=key)
            counter+=1

    plt.xscale('log')
    plt.yscale('log')
    plt.title("{:.2f}GeV Cascade Rates".format(energy/const.GeV))
    plt.xlabel("Parent Neutrino Energy [GeV]")
    plt.ylabel(r"Probability [GeV$^{-1}$]")
    plt.legend()
    print("saving 'wow.png'")
    plt.savefig("wow.png",dpi=400)

savefile = ".analysis_level.dat"
def _load_data():
    """
    Loads the datafile. Returns tuple

    0 - parent energies
    1 - child energies
    3 - muon cascade flux (DD)
    4 - not muon cascade flux (DD)
    """
    print("Loading Data")
    f = open(savefile, 'rb')
    all_data = pickle.load(f)
    f.close()
    return( all_data["parent_energies"], all_data["child_energies"], \
                all_data["muon_ones"], all_data["not_muon"] )
    
def _save_data(parent_energies, child_energies, muon_ones, not_muon):
    """
    Saves the generated data for use later. 
    """
    all_data = {"parent_energies": parent_energies, 
                    "child_energies": child_energies, 
                    "muon_ones": muon_ones, 
                    "not_muon": not_muon}
    f = open(savefile,'wb')
    pickle.dump( all_data, f, -1)
    f.close()   
    print("Data Saved!")

def generate_mode2_data():
    """
    Generates the data for mode 2
    """
    raise DeprecationWarning("This function is deprecated! ")
    print("Calculating E_i Distributions")
    e_min = 50*const.GeV
    e_max = 100*const.TeV

    these_energies = np.logspace(np.log10(e_min), np.log10(e_max), n_bins)  # cascade energy 
    parent_energies = np.logspace(np.log10(e_min), np.log10(e_max)+2,n_bins)
    widths = np.zeros(n_bins)
    cascade_widths = get_width(these_energies)/const.GeV

    #               [cascade energy] [parent energy]
    muon_ones = np.zeros((n_bins, n_bins))
    not_muon  = np.zeros((n_bins, n_bins))

    for index in range(n_bins): #integrate over cascade energies 
        # let's populate those 2D lists
        lepton_energies = parent_energies - these_energies[index] # lepton energy 
        actual_widths = get_width(lepton_energies)/const.GeV
        from_diffy, widths = get_distribs_for_cascade_energy(these_energies[index], actual_widths, parent_energies)
        for flav in flavors:
            for neut in neuts:
                for curr in currents:
                    if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                        continue 
                    key = flav+'_'+neut + '_'+curr
                    if curr=="NC" and flav=="Mu":
                        muon_ones[index] += from_diffy[key]
                    else:
                        not_muon[index]  += from_diffy[key] 
        if do_norm:  
            norm = sum(muon_ones[index]*widths) + sum(not_muon[index]*widths)
        else:
            norm = 1.
        if norm!=0 and norm!=0.0:
            for i in range(n_bins):
                muon_ones[index][i] = muon_ones[index][i]/norm
                not_muon[index][i] = not_muon[index][i]/norm

    _save_data(parent_energies, these_energies, muon_ones, not_muon)
    return( parent_energies, these_energies, muon_ones, not_muon )

# from_diffy[key][j] =data.get_flux(event_energies[j], key)*get_diff_xs(event_energies[j], get_flavor(key), get_neut(key), get_curr(key),lepton_energy[j],0.5)*event_widths[j]
def get_2D_flux(lepton_energies, event_energies, key):
    """
    Builds a doubly differential flux array for the given desired flux (as specified by the key)
    """
    if not (isinstance(lepton_energies, list) or isinstance(lepton_energies, np.ndarray)):
        raise TypeError("Expected {} for 'lepton_energies', got {}".format(np.ndarray, type(lepton_energies)))
    if not (isinstance(event_energies, list) or isinstance(event_energies, np.ndarray)):
        raise TypeError("Expected {} for 'event_energies', got {}".format(np.ndarray, type(event_energies)))
    if not isinstance(key, str):
        raise TypeError("Expected {} for key, got {}".format(str, type(key)))

    lep_bins = len(lepton_energies)
    evt_bins = len(event_energies)
    fluxes = np.zeros((lep_bins,evt_bins))

    for lep_bin in range(lep_bins):
        for evt_bin in range(evt_bins):
            if event_energies[evt_bin]<lepton_energies[lep_bin]:
                continue
            fluxes[lep_bin][evt_bin] = data.get_flux(event_energies[evt_bin],key)
            fluxes[lep_bin][evt_bin]*= get_diff_xs(event_energies[evt_bin], get_flavor(key), get_neut(key), get_curr(key), lepton_energies[lep_bin],0.0) 
    return(fluxes)

def swap_to_cascade(flux, lepton_energies, event_energies, key):
    """
    Takes a doubly differential flux array and the lepton and event energies

    Then it swaps it around to be doubly differential in cascade and event energies 
    """
    # we want this to follow the same trend 
    cascade_energies = np.array([i for i in lepton_energies])
    flav = key.split('_')[0]
    curr = key.split('_')[2]
    cce_eff = 0.98

    # cascade vs event 
    nuflux = np.zeros((len(cascade_energies), len(event_energies)))

    def get_cascade_bin( energy ):
        """
        Returns the appropriate bin for a provided cascade energy 
        """
        min_diff = None
        minbin = None
        for i in range(len(cascade_energies)):
            if min_diff is None:
                min_diff = abs(energy-cascade_energies[i])
                minbin = i
            else:
                if min_diff > abs(energy-cascade_energies[i]):
                    min_diff = abs(energy-cascade_energies[i])
                    minbin = i
                else:
                    return(i)
        return(i)

    lepton_widths = get_width(lepton_energies)/const.GeV
    cascade_widths = get_width(cascade_energies)/const.GeV

    for lep_bin in range(len(lepton_energies)):
        for evt_bin in range(len(event_energies)):
            #if lepton_energies[lep_bin] > event_energies[evt_bin]:
            #    continue #forbidden! 
            if curr=='CC':
                if flav=="E":
                    index = get_cascade_bin(event_energies[evt_bin]*cce_eff)
                else:
                    continue # track
            else: #NC
                index = get_cascade_bin(event_energies[evt_bin] - lepton_energies[lep_bin])
            if index!=len(cascade_energies) and index!=-1:
                nuflux[index][evt_bin] += flux[lep_bin][evt_bin]*lepton_widths[lep_bin]/cascade_widths[index]
            else:
                print("Error? Key {}, event {:.2f} lepton {:.2f}".format(key, event_energies[evt_bin]/const.GeV, lepton_energies[lep_bin]/const.GeV))
    return(flux)

def generate_doubly_diff_fluxes(n_bins=200, debug=False, flavor_split=False):
    """
    Newer method for working out the likely event energy based on the cascade energy 

    Now it prepares a doubly differential flux (lepton vs event)
    Then it transitions to a doubly differential one with respect to the cascade 
        - it uses the lepton and event energy to get the hadronic component (or 98% of the event energy for CC e-)
        - moves the flux from the lepton bin to the corresponding cascade energy bin, and rescales according to the different lepton/cascade bin widths 

    Then it combines these into two 2D differential things for the "muon" ones and "not muon" ones (simple addition)
    and returns the event energies, lepton energies, and the two 2D differential arrays 

    If "debug" is True, it instead returns the og nuflux dictionary 
    """
    print("Generating 2D LeptonVsEvent")
    e_min = 50*const.GeV
    e_max = 100*const.TeV
    extra = 2

    lepton_energies = np.logspace(np.log10(e_min), np.log10(e_max),n_bins)
    event_energies = np.logspace(np.log10(e_min), np.log10(e_max)+extra,n_bins)
    if debug:
        print("lepton ranges from {:.2f} to {:.2f} GeV".format(min(lepton_energies)/const.GeV, max(lepton_energies)/const.GeV))
        print("event ranges from  {:.2f} to {:.2f} GeV".format(min(event_energies)/const.GeV,  max(event_energies)/const.GeV))
    flux = {}
    nuflux = {}
    from_muon = np.zeros((n_bins, n_bins))
    from_not = np.zeros((n_bins, n_bins))
    
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                key = flav+"_"+neut+"_"+curr
                if (flav=="Mu" or flav=="Tau") and curr=="CC":
                    continue # skip tracks
                else:
                    print("Working on {} {} {} Events".format(curr, neut, flav))
                    # [cascade/lepton] [ event ]
                    flux[key] = get_2D_flux( lepton_energies, event_energies, key)
                    nuflux[key] = swap_to_cascade( flux[key], lepton_energies, event_energies,key)
                    if flav=="Mu":
                        from_muon += flux[key] if debug else nuflux[key]
                    else:
                        from_not += flux[key] if debug else nuflux[key]
#parent child muon not_muon
    _save_data( event_energies, lepton_energies, from_muon, from_not)
    if flavor_split:
        return(event_energies,lepton_energies,flavor_split)
    else:
        return(event_energies,lepton_energies, from_muon, from_not) 

if mode==8:
    n_bins = 200
    event_energies, cascade_energies, from_muon, from_not = generate_doubly_diff_fluxes(n_bins)
    
    from_muon = np.ma.masked_where(from_muon<=0, from_muon)
    from_not  = np.ma.masked_where(from_not<=0, from_not)

    plt.figure()
    levels = np.logspace(-70,-45,10)
    print("Max of muon: {}".format(np.min(from_muon)))
    cf = plt.contourf(event_energies/const.GeV, cascade_energies/const.GeV, from_muon,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Parent Energy [GeV]', size=14)
    plt.ylabel('Cascade Energy [GeV]', size=14)
    plt.grid(which='major',alpha=0.7)
    cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
    print("saving from_muon.png")
    plt.savefig('from_muon.png', dpi=400)
    plt.clf()
    cf = plt.contourf(event_energies/const.GeV, cascade_energies/const.GeV, from_not,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
    cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
    plt.grid(which='major',alpha=0.7)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Parent Energy [GeV]', size=14)
    plt.ylabel('Cascade Energy [GeV]', size=14)
    plt.savefig('from_not.png', dpi=400)
    print("Saving from_not.png")


if mode==2 or mode==4:
    """
    Creates two 2D contour plots showing doubly differential event rates as a function of event and cascade energy. 

    Deprecated. See generate_doubly_diff_fluxes! 
    """
    if load_stored and os.path.exists(savefile):
        parent_energies, these_energies, muon_ones, not_muon = _load_data()
    else:
        parent_energies, these_energies, muon_ones, not_muon = generate_doubly_diff_fluxes(n_bins)

    print("Plotting")
    # so it doesn't scream about logged zeros 
    muon_ones = np.ma.masked_where(muon_ones<=0, muon_ones)
    not_muon  = np.ma.masked_where(not_muon<=0, not_muon)
    levels = np.logspace(-5,0,8) if do_norm else np.logspace(-50,-30, 8)
    
    if mode==2:
        # evt /s /cm2 /GeV /sr 
        plt.figure()
        cf = plt.contourf(parent_energies/const.GeV, these_energies/const.GeV, muon_ones,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Parent Energy [GeV]', size=14)
        plt.ylabel('Cascade Energy [GeV]', size=14)
        plt.grid(which='major',alpha=0.7)
        cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
        cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
        print("saving muon_ones.png")
        plt.savefig('muon_ones.png', dpi=400)
        plt.clf()
        cf = plt.contourf(parent_energies/const.GeV, these_energies/const.GeV, not_muon,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
        cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
        cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
        plt.grid(which='major',alpha=0.7)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Parent Energy [GeV]', size=14)
        plt.ylabel('Cascade Energy [GeV]', size=14)
        plt.savefig('not_muon.png', dpi=400)
        print("Saving not_muon.png")
    elif mode==4:
        levels = np.logspace(-2,2,20)
        plt.figure()
        cf = plt.contourf(parent_energies/const.GeV, these_energies/const.GeV, muon_ones/not_muon,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Parent Energy [GeV]', size=14)
        plt.ylabel('Cascade Energy [GeV]', size=14)
        plt.grid(which='major',alpha=0.7)
        cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
        cbar.set_label("Muon Rate / Not Muon Rate")
        print("Saving ratio_plot.png")
        plt.savefig('ratio_plot.png', dpi=400) 
if mode==5:
    """
    This mode prepares a plot showing the median energy of an event as a function of cascade energy.
    It has two trends, one for muons and one for everything else.

    Error bands are shown representing a 1-sigma deviation 
    """
    if load_stored and os.path.exists(savefile):
        parent_energies, these_energies, muon_ones, not_muon = _load_data()
    else:
        parent_energies, these_energies, muon_ones, not_muon = generate_doubly_diff_fluxes(n_bins)
   
    cascade_widths = get_width(these_energies)/const.GeV
    parent_widths = get_width(parent_energies)/const.GeV

    # the [outer_index] refers to an observed cascade energy
    # the [inner_inde] corresponds to (lower_extent_error, mean, upper_extent_error)
    muon_expectation = np.array([(0.,0.,0.) for j in range(n_bins)])
    nute_expectation = np.array([(0.,0.,0.) for j in range(n_bins)])
    p_muon = np.zeros(n_bins)
    for index in range(n_bins): # iterate over the cascade energies 
        scale = 1.# cascade_widths[index]
        mean, sigma_up, sigma_down = get_exp_std( parent_widths, muon_ones[index], parent_energies/const.GeV)
        muon_expectation[index] = np.array([sigma_down, mean, sigma_up])
        
        mean, sigma_up, sigma_down = get_exp_std( parent_widths, not_muon[index], parent_energies/const.GeV)
        nute_expectation[index] = np.array([sigma_down, mean, sigma_up])

        p_muon[index] = sum(muon_ones[index]*parent_widths)/(sum(not_muon[index]*parent_widths) + sum(muon_ones[index]*parent_widths))
        
    figs,axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

    # we need to transpose the expectations for ease of plotting 
    muon_expectation = np.transpose(muon_expectation)
    nute_expectation = np.transpose(nute_expectation)


    axes[0].fill_between( these_energies/const.GeV, muon_expectation[1]-muon_expectation[0], muon_expectation[1]+muon_expectation[2], color='#5f97c7',alpha=0.2)
    axes[0].fill_between( these_energies/const.GeV, nute_expectation[1]-nute_expectation[0], nute_expectation[1]+nute_expectation[2], color='#b31007', alpha=0.2)
    axes[0].plot( these_energies/const.GeV, muon_expectation[1], drawstyle='steps', label="Muon Origin", color='#5f97c7')
    axes[0].plot( these_energies/const.GeV, nute_expectation[1], drawstyle='steps', label="Not Muon", color='#b31007')
    axes[1].plot(these_energies/const.GeV, p_muon)
    
    axes[0].set_xlim([5e1, 10**5])
    axes[1].set_xlim([5e1, 10**5])
    axes[0].set_ylim([10**1, 10**7])
    axes[1].set_ylim([1e-4,1])
    axes[0].grid('major', alpha=0.5 )
    axes[0].legend()

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving probable_energy.png")
    plt.savefig("probable_energy.png", dpi=400)

if mode==6:
    """
    In this mode, I make a two-for-one plot. We show the most probable event energy as a
        function of the cascade's energy. Underneath this plot, we show the probability 
        that the event is due to a muon.
    
    This doesn't make the muon/notmuon distinction from mode 5
    """
    if load_stored and os.path.exists(savefile):
        parent_energies, these_energies, muon_ones, not_muon = _load_data()
    else:
        parent_energies, these_energies, muon_ones, not_muon = generate_doubly_diff_fluxes(n_bins)
  
    cascade_widths = get_width(these_energies)/const.GeV
    parent_widths = get_width(parent_energies)/const.GeV

    expectation = np.array([ np.array([0.,0.,0.]) for j in range(n_bins)])
    p_muon = np.zeros(n_bins)
    for index in range(n_bins):
        norm = sum(muon_ones[index]*parent_widths) + sum(not_muon[index]*parent_widths)
        mean, sigma_up, sigma_down = get_exp_std( parent_widths, muon_ones[index], parent_energies/const.GeV)        
        expectation[index] += np.array([sigma_down, mean, sigma_up])*sum(muon_ones[index]*parent_widths)/norm

        mean, sigma_up, sigma_down = get_exp_std( parent_widths, not_muon[index], parent_energies/const.GeV)
        expectation[index] += np.array([sigma_down, mean, sigma_up])*sum(not_muon[index]*parent_widths)/norm

        p_muon[index] = sum(muon_ones[index]*parent_widths)/norm

    figs,axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

    # we need to transpose the expectations for ease of plotting 
    expectation = np.transpose(expectation)

    print("Performing Linear regression")
    m, b = np.polyfit(x=np.log10(these_energies/const.GeV), y=np.log10(expectation[1]/const.GeV), deg=1)
    print(" E [GeV] = ({:.2f})*(Cascade E/GeV)^{} ".format((1e9)*(10**b), m))

    axes[0].fill_between( these_energies/const.GeV, expectation[1]-expectation[0], expectation[1]+expectation[2], color='#5f97c7',alpha=0.2)
    axes[0].plot( these_energies/const.GeV, expectation[1], drawstyle='steps', color='#5f97c7')
    axes[1].plot(these_energies/const.GeV, p_muon)
    
    axes[0].set_xlim([5e1, 10**5])
    axes[1].set_xlim([5e1, 10**5])
    axes[0].set_ylim([10**1, 10**7])
    axes[1].set_ylim([1e-4,1])
    axes[0].grid('major', alpha=0.5 )

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving predicted_energ_E.png")
    plt.savefig("predicted_event_E.png", dpi=400)

if mode==7:
    """
    This is deprectated.

    It was another method to recreate the results from mode 6, but ultimately encountered the same issues as mode 6.
    """
    e_min = 50*const.GeV
    e_max = 100*const.TeV

    cascade_energies = np.logspace(np.log10(e_min), np.log10(e_max), n_bins)  # cascade energy 
    event_energies = np.logspace(np.log10(e_min), np.log10(e_max)+2,n_bins)
    widths = np.zeros(n_bins)
    cascade_widths = get_width(cascade_energies)/const.GeV


    #               [event energy] [parent energy]

    expectation = np.array([ np.array([0.,0.,0.]) for i in range(n_bins)])
    p_muon = np.zeros(n_bins)

    for index in range(n_bins):
        # let's populate those 2D lists
        muon_ones = np.zeros(n_bins)
        not_muon = np.zeros(n_bins)

        lepton_widths = get_width(event_energies-cascade_energies[index])/const.GeV
        from_diffy, widths = get_distribs_for_cascade_energy(cascade_energies[index],lepton_widths, event_energies)
        for flav in flavors:
            for neut in neuts:
                for curr in currents:
                    if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                        continue 
                    key = flav+'_'+neut + '_'+curr
                    if curr=="NC" and flav=="Mu":
                        muon_ones += from_diffy[key]
                    else:
                        not_muon  += from_diffy[key] 
        norm = sum(muon_ones*widths) + sum(not_muon*widths)
            
        mean, sigma_up, sigma_down = get_exp_std( widths, muon_ones, event_energies/const.GeV)
        expectation[index] += np.array([sigma_down, mean, sigma_up])*sum(muon_ones*widths)/norm
 
        mean, sigma_up, sigma_down = get_exp_std( widths, not_muon, event_energies/const.GeV)
        expectation[index] += np.array([sigma_down, mean, sigma_up])*sum(not_muon*widths)/norm

        p_muon[index] = sum(muon_ones*widths)/norm
    
    figs,axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

    # we need to transpose the expectations for ease of plotting 
    expectation = np.transpose(expectation)

    axes[0].fill_between( cascade_energies/const.GeV, expectation[1]-expectation[0], expectation[1]+expectation[2], color='#5f97c7',alpha=0.2)
    axes[0].plot( cascade_energies/const.GeV, expectation[1], drawstyle='steps', color='#5f97c7')
    axes[1].plot(cascade_energies/const.GeV, p_muon)

    axes[1].grid('major', alpha=0.5)
    axes[0].set_xlim([5e1, 10**5])
    axes[1].set_xlim([5e1, 10**5])
    axes[0].set_ylim([10**1, 10**7])
    axes[1].set_ylim([1e-4,1])
    axes[0].grid('major', alpha=0.5 )

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving predicted_event_E_two.png")
    plt.savefig("predicted_event_E_TWO.png", dpi=400)



if mode==9:
    """
    This mode prepares a plot showing the median energy of an event as a function of cascade energy.
    It has two trends, one for muons and one for everything else.

    Error bands are shown representing a 1-sigma deviation 
    """
    n_bins = 200
    parent_energies, these_energies, flux = generate_doubly_diff_fluxes(n_bins, flavor_split=True)
   
    cascade_widths = get_width(these_energies)/const.GeV
    parent_widths = get_width(parent_energies)/const.GeV

    taus =np.zeros(n_bins)
    eles =np.zeros(n_bins)
    muon =np.zeros(n_bins)

    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                    continue 
                key = flav+'_'+neut + '_'+curr
                if flav=='Mu':
                    muon+=flux[key]
                elif flav=='Tau':
                    taus+=flux[key]
                elif flav=='E':
                    eles+=flux[key]
                else:
                    raise ValueError("Explain yourself. {}".format(flav))

    # the [outer_index] refers to an observed cascade energy
    # the [inner_inde] corresponds to (lower_extent_error, mean, upper_extent_error)
    muon_expectation = np.array([(0.,0.,0.) for j in range(n_bins)])
    taus_expectation = np.array([(0.,0.,0.) for j in range(n_bins)])
    eles_expectation = np.array([(0.,0.,0.) for j in range(n_bins)])

    p_muon = np.zeros(n_bins)
    for index in range(n_bins): # iterate over the cascade energies 
        scale = 1.# cascade_widths[index]
        mean, sigma_up, sigma_down = get_exp_std( parent_widths, muon[index], parent_energies/const.GeV)
        muon_expectation[index] = np.array([sigma_down, mean, sigma_up])
        
        mean, sigma_up, sigma_down = get_exp_std( parent_widths, eles[index], parent_energies/const.GeV)
        eles_expectation[index] = np.array([sigma_down, mean, sigma_up])

        mean, sigma_up, sigma_down = get_exp_std( parent_widths, taus[index], parent_energies/const.GeV)
        taus_expectation[index] = np.array([sigma_down, mean, sigma_up])

        p_muon[index] = sum(muon_ones[index]*parent_widths)/(sum(not_muon[index]*parent_widths) + sum(muon_ones[index]*parent_widths))
        
    figs,axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

    axes[0].fill_between( these_energies/const.GeV, muon_expectation[1]-muon_expectation[0], muon_expectation[1]+muon_expectation[2], color='#5f97c7',alpha=0.2)
    axes[0].fill_between( these_energies/const.GeV, eles_expectation[1]-eles_expectation[0], eles_expectation[1]+eles_expectation[2], color='#b31007', alpha=0.2)
    axes[0].fill_between( these_energies/const.GeV, taus_expectation[1]-taus_expectation[0], taus_expectation[1]+taus_expectation[2], color='#86d1b2', alpha=0.2)
    axes[0].plot( these_energies/const.GeV, muon_expectation[1], drawstyle='steps', label="Muon Origin", color='#5f97c7')
    axes[0].plot( these_energies/const.GeV, eles_expectation[1], drawstyle='steps', label="Electron Origin", color='#b31007')
    axes[0].plot( these_energies/const.GeV, taus_expectation[1], drawstyle='steps', label="Not Muon", color='#86d1b2')
    axes[1].plot(these_energies/const.GeV, p_muon)
    
    axes[0].set_xlim([5e1, 10**5])
    axes[1].set_xlim([5e1, 10**5])
    axes[0].set_ylim([10**1, 10**7])
    axes[1].set_ylim([1e-4,1])
    axes[0].grid('major', alpha=0.5 )
    axes[0].legend()

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    plt.savefig("all_three_flavor.png",dpi=400)
    print("Saving all_three_flavor.png")
    
