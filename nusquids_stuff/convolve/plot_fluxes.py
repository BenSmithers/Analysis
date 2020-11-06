#!/usr/bin/python3.6
'''
This script plots the fluxes output by the convolve cpp script. It also does some analysis and saves some arrays to disk 

But how does it work?
 1. Raw flux data from the included mceq+nuSQuIDS flux is loaded in by the Data Object. This object has functionality for sampling from the flux at arbitrary parent neutrino energy by interpolating between nsq points.
 2. This flux data is convolved with some cross sections to make a singly-differential flux array of deposited energy vs event energy 

The heavy lifting goes on in the generate_singly_diff_fluxes function 

Different "modes" are used to make different plots. 

Notes:
    + energies are generally in eV for consistency with nuSQuIDs
    + widths are in GeV for consistency with nuSQuIDS (why...)
    + the fluxes output by the Data object are in units of [GeV cm2 s]^-1 (incoming neutrino energy)
    + the cross sections are in units of cm2/GeV (outgoing neutrino energy)
'''

from optparse import OptionParser
import sys
mode_str = "\nModes\n\
Debug - Raw fluxes for the three flavors. No cross sections\n \
4 - PLOT 2d HIST OF parent vs cascade, Ratio plot\n \
5 - The plot showing most probable energy for muon events and non-muon events\n \
6 - The plot just showing must probable energy for all of them\n \
8 - Makes two 2D histograms of parent vs cascade fluxes\n \
9 - Plot of Median event energy for all three flavors\n \
"

parser = OptionParser()
#mode_str = "0 - plot muon vs not muon\n1 - plot all the keys\n 2 - plot 2D hist of parent vs cascade"
parser.add_option("-m", "--mode",
                dest="mode",
                default="0",
                type=str,
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
parser.add_option("-n", "--nbins",
                dest="n_bins",
                default=200,
                type=int,
                help="Number of bins to use for each axis")
parser.add_option("-a", "--angle",
                dest="angle",
                default=None,
                type=float,
                help="At which angle should the plots be made?")

options, args = parser.parse_args()
mode = options.mode
load_stored = options.load_stored
debug = options.debug
glob_angle = options.angle

n_bins = options.n_bins
do_norm=False # deprecated 

if mode.lower()=='a' or mode.lower()=='all':
    mode = 8
    do_all = True
    load_stored = True
else:
    do_all = False
    mode = int(mode)
    recognized_modes = [0,1,2,3,4,5,6,7,-1,8,9]
    if mode not in recognized_modes: 
        raise ValueError("Unrecognized Mode: {}".format(mode))
    if mode in [-1,1,2,3,7]:
        raise DeprecationWarning("Mode {} is deprecated".format(mode))

print("Configuration...")
print("    In Mode {}".format(mode))
print("    Will Load Data" if load_stored else "    Will Generate Data")
print("    Using {} bins".format(n_bins))

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
from utils import get_flavor, get_neut, get_curr, get_exp_std, get_width, get_nearest_entry_to
from utils import bhist
from utils import Data, get_index, get_loc

# tau stuff
from tau_funcs import TauData

# reconstruction data
from deporeco import DataReco

const = nsq.Const()
# colormap thing
cmap = plt.get_cmap('coolwarm')
n_colors = 6
def get_color(which, how_many=n_colors):
    return( cmap( float(which)/how_many ) )

# load the data using the default filename, 'atmosphere.txt'
data = Data()
tauData = TauData()

scale_e = np.array(data.energies)

# this block here is supposed to just plot all the raw fluxes 
if debug:
    print("Making Debug plot")
    binny = len(scale_e)
    taus = np.zeros(binny)
    muon = np.zeros(binny)
    ele = np.zeros(binny)

    angle_bin = 0
    print("Angle: {}".format(data.angles[angle_bin]))

    for key in data.fluxes:
        flav=str(key).split('_')[0]
        if 'Tau'==flav:
            taus+=[data.fluxes[key][i][angle_bin] for i in range(len(data.fluxes[key]))]
        elif "E"==flav:
            ele+=[data.fluxes[key][i][angle_bin] for i in range(len(data.fluxes[key]))]
        elif "Mu"==flav:
            muon+=[data.fluxes[key][i][angle_bin] for i in range(len(data.fluxes[key]))]
        else:
            raise Exception("You might have done steriles? {} is unrecognized".format(flav))


    plt.plot( scale_e/const.GeV, muon, color=get_color(0,3),label="Muons")
    plt.plot(scale_e/const.GeV, taus, color=get_color(1,3),label="Taus")
    plt.plot(scale_e/const.GeV, ele, color=get_color(2,3), label="Electrons")
    plt.legend()
    plt.xlabel("Event Energy [GeV]",size=14)
    plt.ylabel("Flux [GeV sec cm$^{2}$ sr]$^{-1}$",size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('raw_fluxes_{:.2f}.png'.format(glob_angle), dpi=400)
    print("saved raw fluxes")
    plt.close()


if mode==1:
    """
    This was used to make plots showing the probability density of event energy using a fixed cascade energy.

    This is deprecated. 
    """

    e_min = 10*const.GeV
    e_max = 10*const.PeV
    in_energies = np.logspace(np.log10(e_min), np.log10(e_max), n_bins)
    # the initial neutrino energy we are considering

    energy = 100*const.GeV
    casc_widths = get_width(in_energies - energy)
    from_diffy, widths  = get_distribs_for_cascade_energy(energy, casc_widths/const.GeV, in_energies)
    from_muon = np.array([ 0. for ienergy in range(n_bins) ])
    from_not = np.array([ 0. for ienergy in range(n_bins) ])
    for flav in data.flavors:
        for neut in data.neuts:
            for curr in data.currents:
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
    plt.savefig("wow_{:.2f}.png".format(glob_angle),dpi=400)

savefile = ".analysis_level.dat"
def _load_data():
    """
    Loads the datafile. Returns tuple

    0 - parent energies
    1 - child energies
    2 - nuflux
    3 - cos(zenith) edges
    """
    print("Loading Data")
    f = open(savefile, 'rb')
    all_data = pickle.load(f)
    f.close()
    angle_edges = all_data["angle_edges"]
    nuflux = all_data["nuflux"]

    if glob_angle is not None:
        # figure out which slice to return
        # [:,:,N] would return the N'th angle bin
        if glob_angle<min(angle_edges) or glob_angle>max(angle_edges):
            for key in nuflux.keys():
                nuflux[key] = nuflux[key][:,:,0]*0.
        else:
            lower, upper = get_loc(glob_angle, angle_edges)

            print("Grabbed angle bin {}".format(lower))
            width = abs( np.arccos(angle_edges[upper]) - np.arccos(angle_edges[lower]))

            for key in nuflux.keys():
                nuflux[key] = nuflux[key][:,:,lower]*width
                

    return( all_data["parent_energies"], all_data["child_energies"], \
                nuflux, angle_edges )
    
def _save_data(e_reco, e_true, a_reco, a_true, flux):
    """
    Saves the generated data for use later. 
    """
    all_data = {"e_reco": e_reco,
                "e_true": e_true,
                "a_reco": a_reco,
                "a_true": a_true,
                "flux": flux}
    f = open(savefile,'wb')
    pickle.dump( all_data, f, -1)
    f.close()   
    print("Data Saved!")



def do_for_key(event_edges,cascade_edges, key, angles=None):
    """
    This function takes the desired bin edges for the event energies and deposited energies along with the dictionary key corresponding to a specific combination of falvor, current, and neutrino type.

    It builds up the 2D flux array (singly differential), which it then returns 
    """
    evt = bhist([event_edges])
    cas = bhist([cascade_edges])
    reco = bhist([cascade_edges])

    event_energies = evt.centers
    event_widths = evt.widths
    cascade_energies = cas.centers
    cascade_widths = cas.widths
    reco_energies = reco.centers
    reco_widths = reco.widths 

    flav = key.split("_")[0]
    curr = key.split("_")[2]

    if angles is None:
        ang_list = [None]
    else:
        ang_list = angles

    if angles is None:
        flux = bhist((cascade_edges, event_edges))
    else:
        flux = bhist((cascade_edges, event_edges, angles))

    for angle in ang_list:
        # need special Tau treatment 
        if curr=="CC":
            # deposit all the energy. Hadronic and Leptonic (event) contribute 
            # Since the energy always all gets deposited, we only need to do one loop!
            # So, for a given "deposited energy" (cascade_energy), we already know the total energy. 
            # Therefore we just get the total cross section * flux there... the units end up as [s GeV in sr]^-1 
            for cas_bin in range(len(cascade_energies)): 
                deposited_energy = cascade_energies[cas_bin] # energy going to the leptonic component! 
                 
                if flav.lower()=='tau':
                    # Etau is cascade_energies[cas_bin]
                    # How much energy is visible in the various tau decays though? 
                    # going from zero->deposited energy
                    deposited_energy = 0.5*(tauData(deposited_energy/const.GeV,1)+tauData(deposited_energy/const.GeV,-1))

                amount =data.get_flux(deposited_energy,key, angle=angle)
                amount *= get_diff_xs(deposited_energy, get_flavor(key), get_neut(key), get_curr(key))
                try:
                    if angle is None:
                        flux.register( amount, cascade_energies[cas_bin], deposited_energy) # add it in! 
                    else: 
                        flux.register( amount, cascade_energies[cas_bin], deposited_energy, angle)
                except ValueError:
                    if flav.lower()!='tau':
                        raise Exception("It wasn't tau. Something's wrong")

        else:
            # in this case, knowing the cascade doesn't tell us anything about the event energy. 
            # so we loop over both, get the flux*differential_xs at each bin combination, and multiply by the widths of deposited-energy-bin to get the same units as in the CC case 
            for evt_bin in range(len(event_energies)):
                for cas_bin in range(len(cascade_energies)):
                    lepton_energy = event_energies[evt_bin] - cascade_energies[cas_bin]

                    # we'll have nowhere to put these, so let's just skip this
                    if lepton_energy < min(cascade_energies):
                        continue
                    if lepton_energy > max(cascade_energies):
                        continue

                    amount =data.get_flux(event_energies[evt_bin],key, angle=angle)
                    amount *= get_diff_xs(event_energies[evt_bin], get_flavor(key), get_neut(key), get_curr(key), lepton_energy,0.0)*cascade_widths[cas_bin]
                    if angle is None:
                        flux.register(amount, cascade_energies[cas_bin], event_energies[evt_bin])
                    else:
                        flux.register(amount, cascade_energies[cas_bin], event_energies[evt_bin], angle)
    
    # build a new bhist in reconstruction space (Still with event energy too)
    # then scan through deposited-true angle space
    # and use the PDFS to smear the true values into reconstructed values, depositing them into the reconstruction bhist  

    return(flux.fill)

def generate_singly_diff_fluxes(n_bins,debug=False):
    """
    This is the newest, most streamlined function for generating the singly-differential flux arrays. 
    It has the same functionality as the old ones, but now it skips the step of filling in the regular flux array before swapping to the deposited energy flux array. 

    It goes over the different kind of fluxes and builds up a 2D or 3D array, convolving the entries with the relevant differential cross section
        the dimensionality depends on whether or we are integrating over the zenith angles
    """
    e_min = 10*const.GeV
    e_max = 100*const.TeV
    extra = 2
    
    all_angles = data.angles
    
    event_edges = np.logspace(np.log10(e_min), np.log10(e_max)+extra,n_bins+1)
    cascade_edges = np.logspace(np.log10(e_min), np.log10(e_max),n_bins+1)
    angle_edges = np.linspace(min(all_angles), max(all_angles), n_bins+1) # angle is in cos(zenith), so like -1->0

    sep_angles = True
    if sep_angles:
        from_muon = np.zeros((n_bins,n_bins,n_bins))
        from_not = np.zeros((n_bins,n_bins,n_bins))
    else:
        from_muon = np.zeros((n_bins,n_bins))
        from_not = np.zeros((n_bins,n_bins))

    nuflux = {}

    for key in data.get_keys(): #flavor, current, interaction 
        nuflux[key] = do_for_key(event_edges,cascade_edges,key, (angle_edges if sep_angles else None))

    # if global variable "angle" isn't none, then we can separate out just a single angle
    
    if False: # glob_angle is not None:
        # figure out which slice to return
        # [:,:,N] would return the N'th angle bin
        if glob_angle<min(angle_edges) or glob_angle>max(angle_edges):
            for key in nuflux.keys():
                nuflux[key] = nuflux[key][:,:,0]*0.
        else:
            # get the appropriate entry for this angle and extract the slice 
            lower, upper = get_loc(glob_angle, angle_edges)
            for key in nuflux.keys():
                nuflux[key] = nuflux[key][:,:,lower]

    return(event_edges,cascade_edges, nuflux, angle_edges)

def incorporate_recon(event_edges, cascade_edges, nuflux, angle_edges):
    e_min = min(cascade_edges)
    e_max = max(cascade_edges)

    z_min = min(angle_edges)
    z_max = max(angle_edges)

    cascade_centers = bhist([cascade_edges]).centers
    true_e_centers = bhist([event_edges]).centers
    true_ang_centers = bhist([angle_edges]).centers
    
    r_energy = bhist([ np.logspace(np.log10(e_min), np.log10(e_max), len(cascade_edges)) ])
    r_angle  = bhist([ np.linspace( z_min, z_max, len(angle_edges))])
    
    r_energy_centers = r_energy.centers
    r_energy_widths = r_energy.widths
    r_angle_centers = r_angle.centers
    r_angle_widths = r_angle.widths

    #build the data object
    dataobj = DataReco(r_energy.edges, r_angle.edges, cascade_edges, angle_edges)

    # may god have mercy on our souls 
    recoflux = {}
    for key in nuflux.keys():
        print("Reconstructing {} Flux".format(key))
        # energy x, angle y
        recoflux[key] = np.zeros((len(r_energy_centers),len(true_e_centers), len(r_angle_centers),len(true_ang_centers)))
        for i_e_reco in range(len(r_energy_centers)):
            for i_e_depo in range(len(cascade_centers)):
                depo_odds = dataobj.get_energy_reco_odds(i_e_depo, i_e_reco) #per 
                for i_a_true in range(len(true_ang_centers)):
                    for i_a_reco in range(len(r_angle_centers)):
                        ang_odds = dataobj.get_czenith_reco_odds(i_a_true, i_a_reco) #per sr
                        for i_e_true in range(len(true_e_centers)):
                            amt = nuflux[key][i_e_true][i_e_depo][i_a_true]*depo_odds*ang_odds #per angle per gev depo
                            recoflux[key][i_e_reco][i_e_true][i_a_reco][i_a_true] += amt

    _save_data(r_energy.edges, event_edges, angle_edges,angle_edges, recoflux)
    


def sep_by_flavor(nuflux):
    """
    So this takes that nuflux object, a dictionary of 3D arrays, and separates it into two 3D arrays: one for muons and one for non-muons
    """

    if not isinstance(nuflux, dict):
        raise TypeError("nuflux should be a {}, not a {}".format(dict, type(nuflux)))
    if not isinstance(nuflux[list(nuflux.keys())[0]], np.ndarray):
        raise TypeError("Entries in nuflux should all be {}, not {}".format(np.ndarray, type(nuflux[list(nuflux.keys())[0]])))

    entry_shape = np.shape(nuflux[list(nuflux.keys())[0]])
    from_muons = np.zeros(shape=entry_shape)
    from_not = np.zeros(shape=entry_shape)

    for key in nuflux:
        flavor = key.split('_')[0]
        if flavor=="Mu":
            from_muons+=nuflux[key]
        else:
            from_not+=nuflux[key]
    return(from_muons, from_not)

if (do_all and not load_stored) or mode==0:
    a,b,c,d = generate_singly_diff_fluxes(n_bins)
    incorporate_recon(a,b,c,d)

if mode==8 or do_all:
    if load_stored and os.path.exists(savefile):
        event, cascade, nuflux, angle_edges = _load_data()
    else:
        event, cascade, nuflux, angle_edges = generate_singly_diff_fluxes(n_bins)

    from_muon, from_not = sep_by_flavor(nuflux)

    event_energies = np.array(bhist([event]).centers)
    cascade_energies = np.array(bhist([cascade]).centers)

    from_muon = np.ma.masked_where(from_muon<=0, from_muon)
    from_not  = np.ma.masked_where(from_not<=0, from_not)

    plt.figure()
    levels = np.logspace(-50,-33,10)
    print("Max of muon: {}".format(np.max(from_muon)))
    cf = plt.contourf(event_energies/const.GeV, cascade_energies/const.GeV, from_muon,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Parent Energy [GeV]', size=14)
    plt.ylabel('Cascade Energy [GeV]', size=14)
    plt.grid(which='major',alpha=0.7)
    cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
    print("saving from_muon.png")
    plt.savefig('from_muon_{:.2f}.png'.format(glob_angle), dpi=400)
    plt.clf()
    cf = plt.contourf(event_energies/const.GeV, cascade_energies/const.GeV, from_not,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
    cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
    plt.grid(which='major',alpha=0.7)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Parent Energy [GeV]', size=14)
    plt.ylabel('Cascade Energy [GeV]', size=14)
    plt.savefig('from_not_{:.2f}.png'.format(glob_angle), dpi=400)
    print("Saving from_not.png")


if mode==2 or mode==4 or do_all:
    """
    Creates two 2D contour plots showing doubly differential event rates as a function of event and cascade energy. 

    I think this is deprecated now?
    """
    if load_stored and os.path.exists(savefile):
        parent, these, nuflux, angle_edges  = _load_data()
    else:
        parent, these, nuflux, angle_edges  = generate_singly_diff_fluxes(n_bins)
    
    
    muon_ones, not_muon = sep_by_flavor(nuflux)
    
    parent_energies = np.array(bhist([parent]).centers)
    these_energies = np.array(bhist([these]).centers)

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
        plt.savefig('muon_ones_{:.2f}.png'.format(glob_angle), dpi=400)
        plt.clf()
        cf = plt.contourf(parent_energies/const.GeV, these_energies/const.GeV, not_muon,cmap=cm.coolwarm, locator=ticker.LogLocator(), levels=levels)
        cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
        cbar.set_label(r"$dN/(dE_{f}dE_{i})$ [s$^{-1}$GeV$^{-2}$]")
        plt.grid(which='major',alpha=0.7)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Parent Energy [GeV]', size=14)
        plt.ylabel('Cascade Energy [GeV]', size=14)
        plt.savefig('not_muon_{:.2f}.png'.format(glob_angle), dpi=400)
        print("Saving not_muon.png")
    elif mode==4 or do_all:
        levels = np.logspace(-2,2,11)
        plt.figure()
        cf = plt.contourf(parent_energies/const.GeV, these_energies/const.GeV, muon_ones/not_muon,cmap=cm.coolwarm, locator=ticker.LogLocator(),levels=levels)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Parent Energy [GeV]', size=14)
        plt.ylabel('Cascade Energy [GeV]', size=14)
        plt.grid(which='major',alpha=0.7)
        cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
        cbar.set_label("Muon Rate / Not Muon Rate")
        print("Saving ratio_plot.png")
        plt.savefig('ratio_plot_{:.2f}.png'.format(glob_angle), dpi=400) 
if mode==5 or do_all:
    """
    This mode prepares a plot showing the median energy of an event as a function of cascade energy.
    It has two trends, one for muons and one for everything else.

    Error bands are shown representing a 1-sigma deviation 
    """
    if load_stored and os.path.exists(savefile):
        parent, these, nuflux ,angle_edges = _load_data()
    else:
        parent, these, nuflux, angle_edges = generate_singly_diff_fluxes(n_bins)

    muon_ones, not_muon = sep_by_flavor(nuflux)


    parent_con = bhist([parent])
    these_con = bhist([these])

    cascade_widths = np.array( these_con.widths )/const.GeV
    parent_widths = np.array(parent_con.widths)/const.GeV

    parent_energies=np.array(parent_con.centers)
    these_energies = np.array(these_con.centers)

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
    axes[1].set_ylim([0.5,1])
    axes[1].yaxis.set_ticks(np.linspace(0,1,6))
    axes[0].grid('major', alpha=0.5 )
    axes[0].legend()

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving probable_energy.png")
    plt.savefig("probable_energy_{:.2f}.png".format(glob_angle), dpi=400)

if mode==6 or do_all:
    """
    In this mode, I make a two-for-one plot. We show the most probable event energy as a
        function of the cascade's energy. Underneath this plot, we show the probability 
        that the event is due to a muon.
    
    This doesn't make the muon/notmuon distinction from mode 5
    """
    if load_stored and os.path.exists(savefile):
        parent, these, nuflux, angle_edges = _load_data()
    else:
        parent, these, nuflux, angle_edges = generate_singly_diff_fluxes(n_bins)

    muon_ones, not_muon = sep_by_flavor(nuflux)

    parent_con = bhist([parent])
    these_con = bhist([these])

    parent_energies = np.array(parent_con.centers)
    these_energies = np.array(these_con.centers)

    cascade_widths = np.array(these_con.widths)/const.GeV
    parent_widths = np.array(parent_con.widths)/const.GeV

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
    axes[1].set_ylim([0.5,1])
    axes[1].yaxis.set_ticks(np.linspace(0,1,6))
    axes[0].grid('major', alpha=0.5 )

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    print("saving predicted_energ_E.png")
    plt.savefig("predicted_event_E_{:.2f}.png".format(glob_angle), dpi=400)


if mode==9 or do_all:
    """
    This mode prepares a plot showing the median energy of an event as a function of cascade energy.
    It has separate trends for the elecrons, muons, and taus 

    Error bands are shown representing a 1-sigma deviation 
    """
    
    if load_stored and os.path.exists(savefile):
        parent, these, flux, angle_edges = _load_data()
    else:
        parent, these, flux, angle_edges = generate_singly_diff_fluxes(n_bins)
   
    these_con = bhist([these])
    parent_con = bhist([these])

    parent_energies = np.array(parent_con.centers)
    these_energies = np.array(these_con.centers)

    cascade_widths = np.array(these_con.widths)/const.GeV
    parent_widths = np.array(parent_con.widths)/const.GeV

    taus =np.zeros((n_bins,n_bins))
    eles =np.zeros((n_bins,n_bins))
    muon =np.zeros((n_bins,n_bins))

    for flav in data.flavors:
        for neut in data.neuts:
            for curr in data.currents:
                if flav=='Mu' and curr=='CC': # skip the tracks 
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

        p_muon[index] = sum(muon[index]*parent_widths)/(sum(taus[index]*parent_widths) + sum(muon[index]*parent_widths) + sum(eles[index]*parent_widths))
       
    muon_expectation = muon_expectation.transpose()
    eles_expectation = eles_expectation.transpose()
    taus_expectation = taus_expectation.transpose()

    figs,axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

    axes[0].fill_between( these_energies/const.GeV, muon_expectation[1]-muon_expectation[0], muon_expectation[1]+muon_expectation[2], color='#5f97c7',alpha=0.2)
    axes[0].fill_between( these_energies/const.GeV, eles_expectation[1]-eles_expectation[0], eles_expectation[1]+eles_expectation[2], color='#b31007', alpha=0.2)
    axes[0].fill_between( these_energies/const.GeV, taus_expectation[1]-taus_expectation[0], taus_expectation[1]+taus_expectation[2], color='#08a605', alpha=0.2)
    axes[0].plot( these_energies/const.GeV, muon_expectation[1], drawstyle='steps', label="Muon", color='#5f97c7')
    axes[0].plot( these_energies/const.GeV, eles_expectation[1], drawstyle='steps', label="Electron", color='#b31007')
    axes[0].plot( these_energies/const.GeV, taus_expectation[1], drawstyle='steps', label="Tau", color='#08a605')
    axes[1].plot(these_energies/const.GeV, p_muon)
    
    axes[0].set_xlim([5e1, 10**5])
    axes[1].set_xlim([5e1, 10**5])
    axes[0].set_ylim([10**1, 10**7])
    axes[1].set_ylim([0.5,1])
    axes[1].yaxis.set_ticks(np.linspace(0,1,6))
    axes[0].grid('major', alpha=0.5 )
    axes[0].legend()

    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')

    axes[1].set_xlabel("Cascade Energy [GeV]")
    axes[0].set_ylabel("Median Event Energy [GeV]")
    axes[1].set_ylabel("Probability Muon")
    plt.savefig("all_three_flavor_{:.2f}.png".format(glob_angle),dpi=400)
    print("Saving all_three_flavor.png")
    
