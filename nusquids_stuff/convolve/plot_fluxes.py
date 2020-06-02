#!/usr/bin/python3.6
'''
This script plots the fluxes output by the convolve cpp script

But how does it work?
 1. Raw flux data from the included mceq+nuSQuIDS flux is loaded in by the Data Object. This object has functionality for sampling from the flux at arbitrary parent neutrino energy by interpolating between nsq points.
 2. get_distribs_for_initial_energy is then used to get the spectra of events that could produce a cascade of a specified energy. This uses nuSQuIDS DIS cross sections
 3. These two are used in tandem to make most plots. The generate_mode2_data function is used to make a 2D doubly differential flux array. 

Different "modes" are used to make different plots. 
The "Norm" arg is used to switch between normalizing the flux and and showing the raw Flux rates
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
parser.add_option("-n", "--norm",
                dest="norm",
                default=False,
                action="store_true",
                help="Should it be normed? Flag for True")
parser.add_option("-l", "--load_stored",
                dest="load_stored",
                default=False,
                action="store_true",
                help="Should I try to load stored data rather than regenerate it?")
options, args = parser.parse_args()
mode = options.mode
do_norm = options.norm
load_stored = options.load_stored

recognized_modes = [0,1,2,3,4,5,6,7]
if mode not in recognized_modes:
    raise ValueError("Unrecognized Mode: {}".format(mode))

debug = False

if mode==3:
    mode = 2
    do_norm = False
if mode==5:
    do_norm = True

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
    flavor = split[0]
    variety = split[1] # nu or nu-bar 

    flav_index = flavors.index(flavor)
    variety_index = neuts.index(variety)
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

        # this funnny indexing is a result of the way I output the data from nuSQuIDS
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

        # logarithmic interpolation
        # linear in log-energy space 
        y2 = self.fluxes[key][upper_boundary]
        y1 = self.fluxes[key][lower_boundary]
        x2 = self.energies[upper_boundary]
        x1 = self.energies[lower_boundary]
        
        flux_value = energy*((y2-y1)/(x2-x1) ) + y2 -x2*((y2-y1)/(x2-x1))
        return(flux_value)

data = Data()

scale_e = np.array(data.energies)

# this block here is supposed to just plot all the raw fluxes 
if debug:
    it = 0
    for key in data.fluxes:
        plt.plot( scale_e, data.fluxes[key], color=get_color(it),label=key )
        it+=1

    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('reallycool.png', dpi=400)
    plt.clf()
    
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

n_bins = 200
def get_distribs_for_cascade_energy(event, event_widths, in_energies):
    """
    Returns an un-normalized dictionary of fluxes for each of the progenitor particles.
    The thing returns in units of /GeV/s (energy of outgoing particle)

    Param "event" refers to the measured energy of the cascade in eV. FLOAT 
    Param "event_width" refers to the width of the bin for the event's energy. Needed to properly scale rates. [GEV]
    Param "in_energies" is the list of possible progenitor neutrino energies to consider. in eV. LIST
    """
    if not (isinstance(event, float) or isinstance(event,int)):
        raise TypeError("Expected number-like, not {}".format(type(event)))
    if not (isinstance(in_energies, list) or isinstance(in_energies, np.ndarray) or isinstance(in_energies, tuple)):
        raise TypeError("Expected list-like, not {}".format(type(in_energies)))
    rising = in_energies[1] > in_energies[0]

    CC_E_energy = 0.98 # ratio of energy recovered in a CC electron cascade
    # the initial neutrino energy we are considering  -- commented since this is part of the func args now
    widths = get_width(in_energies)/const.GeV
    #debug
    #print("In Energy Increasing" if in_energies[1]>in_energies[0] else "In Energy Decreasing")
    from_diffy = {}
    for flav in flavors:
        for neut in neuts:
            for curr in currents:
                if (flav=='Mu' or flav=='Tau') and curr=='CC': # skip the tracks 
                    continue                    
                # in-energy needs to be changed. 
                # if it's a normal hadronic cascade (muon NC)
                key = flav+'_'+neut + '_'+curr
                from_diffy[key] = np.array([0.0 for j in range(n_bins)])
                if flav=='E' and curr=='CC':
                    # in this case, we recover the energy from the entire cascade, so we essentially exactly know the energy of the parent neutrino
                    # and so get the total cross-section/flux from that energy and put it all in one bin
                    force_initial_energy = event/CC_E_energy 
                    which_bin = 0 if rising else n_bins-1
                    while in_energies[which_bin]<force_initial_energy:
                        which_bin = which_bin+1 if rising else which_bin-1
                        if (which_bin==n_bins or which_bin==-1):
                            break                
                    if (which_bin!=n_bins and which_bin>=0):
                        # in this case, there is only one parent particle energy which could cause this cascade. So we get the total cross section!
                        # N/(cm2 s GeV )     cm2
                        from_diffy[key][which_bin] = data.get_flux(in_energies[which_bin], key)*get_diff_xs(in_energies[which_bin], get_flavor(key), get_neut(key), get_curr(key))
                else:
                    # N/(cm2 s GeV)     cm2/GeV  
                    for j in range(n_bins):
                        if in_energies[j]-event>0:
                            from_diffy[key][j] =data.get_flux(in_energies[j], key)*get_diff_xs(in_energies[j], get_flavor(key), get_neut(key), get_curr(key),in_energies[j]-event,0.)*abs(event_widths[j])
    return(from_diffy, widths)

if mode==0 or mode==1:

    e_min = 50*const.GeV
    e_max = 10*const.PeV
    in_energies = np.logspace(np.log10(e_min), np.log10(e_max), n_bins)
    # the initial neutrino energy we are considering

    energy = 100*const.GeV
    energy_width = 1*const.MeV
    from_diffy, widths  = get_distribs_for_cascade_energy(energy, energy_width/const.GeV, in_energies)
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

        plt.plot(in_energies/const.GeV, from_muon, color=get_color(0,2),label="Muon Origin")
        plt.plot(in_energies/const.GeV, from_not, color=get_color(1,2),label="Other Origin")
        
        print("Total Muon probability: {}".format( sum( from_muon*widths )))
        print("Total Not Muon probability: {}".format(sum(from_not*widths)))

    elif mode==1:
        norm = sum([ sum(from_diffy[key]) for key in from_diffy ])
        n_colors = len(list(from_diffy.keys()))
        counter = 0
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

def generate_mode2_data():
    """
    Generates the data for mode 2
    """
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

if mode==2 or mode==4:
    if load_stored and os.path.exists(savefile):
        parent_energies, these_energies, muon_ones, not_muon = _load_data()
    else:
        parent_energies, these_energies, muon_ones, not_muon = generate_mode2_data()

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
    if load_stored and os.path.exists(savefile):
        parent_energies, these_energies, muon_ones, not_muon = _load_data()
    else:
        parent_energies, these_energies, muon_ones, not_muon = generate_mode2_data()
   
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
    axes[0].set_ylabel("Likely Parent Neutrino Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving probable_energy.png")
    plt.savefig("probable_energy.png", dpi=400)

if mode==6:
    """
    In this mode, I make a two-for-one plot. We show the most probable event energy as a
        function of the cascade's energy. Underneath this plot, we show the probability 
        that the event is due to a muon.

    There's also a bug in here... 
    """
    if load_stored and os.path.exists(savefile):
        parent_energies, these_energies, muon_ones, not_muon = _load_data()
    else:
        parent_energies, these_energies, muon_ones, not_muon = generate_mode2_data()
  
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

    axes[0].fill_between( these_energies/const.GeV, expectation[1]-expectation[0], expectation[1]+expectation[2], color='#5f97c7',alpha=0.2)
    axes[0].plot( these_energies/const.GeV, expectation[1], drawstyle='steps', label="Predicted Energy", color='#5f97c7')
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
    axes[0].set_ylabel("Likely Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving predicted_energ_E.png")
    plt.savefig("predicted_event_E.png", dpi=400)

if mode==7:
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
    axes[0].plot( cascade_energies/const.GeV, expectation[1], drawstyle='steps', label="Predicted Energy", color='#5f97c7')
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
    axes[0].set_ylabel("Likely Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving predicted_event_E_two.png")
    plt.savefig("predicted_event_E_TWO.png", dpi=400)


