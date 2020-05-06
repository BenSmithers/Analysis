#!/usr/bin/python3.6m

# Ben Smithers
# benjamin.smithers@mavs.uta.edu

# Makes plots from the output of my analysis script. 
# separates things specified with the separate flags but puts them in the same figure

import h5py

from mpl_toolkits.mplot3d import Axes3D
# used for parsing inputs
import sys
print("Running Python version "+sys.version)
from optparse import OptionParser
# used for verifying files exist
import os
import numpy as np
import matplotlib
print("Running Matplotlib version: "+matplotlib.__version__)
## matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import warnings

import json # loads the IceCube Data file 
plt.style.use('/home/benito/Desktop/testing/paper.mplstyle')

# this is required for running as a condor job
os.environ['HDF5_USE_FILE_LOCKING'] ='FALSE'

# disable weighting 
no_weight = False

pi = np.pi

parser = OptionParser()
parser.add_option("-i", "--ops", 
                  dest="ops_file", 
                  default="", 
                  type="str",
                  help="Folder of files to be processed")
parser.add_option("-c","--colorbar",
                  dest="color",
                  default="spring",
                  type="str",
                  help="Which colorbar should I plot with?")
parser.add_option("-x", "--bjorken_x", 
                  dest="do_bjorken_x",
                  action="store_false",
                  default=True,
                  help="Should I skip Bjorken_X?")
parser.add_option("-y", "--bjorken_y", 
                  dest="do_bjorken_y",
                  action="store_false",
                  default=True,
                  help="Should I skip Bjorken_Y?")
parser.add_option("-a", "--hAdrons", 
                  dest="hadrons_cosz",
                  action="store_false",
                  default=True,
                  help="Should I skip Hadrons cosz?")
parser.add_option("-o", "--opening", 
                  dest="opening_angle",
                  action="store_false",
                  default=True,
                  help="Should I skip the opening angle?")
parser.add_option("-s", "--secondary_E", 
                  dest="secondary_e",
                  action="store_false",
                  default=True,
                  help="Should I skip the secondary energy spectrum?")
parser.add_option("-d","--hadrons_def",
                  dest="hadrons_def",
                  action="store_false",
                  default=True,
                  help="Should I skip the hadron deflection?")
parser.add_option("-e","--energy",
                  dest="energy",
                  action="store_false",
                  default=True,
                  help="Should I skip the energy?")
parser.add_option("-b","--BjYvE",
                  dest="BjYvE",
                  action="store_false",
                  default=True,
                  help="Should I skip Bjorken Y vs E?")
parser.add_option("--sep_flavor",
                  dest="sep_flavor",
                  action="store_true",
                  default=False,
                  help="Separate plots by flavor?")
parser.add_option("--sep_matter",
                  dest="sep_matter",
                  action="store_true",
                  default=False,
                  help="Separate plots by flavor?")
parser.add_option("--sep_current",
                  dest="sep_current",
                  action="store_true",
                  default=False,
                  help="Separate plots by flavor?")

options, args = parser.parse_args()

# using lots of args to wrest control over how the plots are made
#store parser args! 
sep_flavor    = options.sep_flavor
sep_current   = options.sep_current
sep_matter    = options.sep_matter
ops_file    = options.ops_file
color         = options.color
bjorken_x     = options.do_bjorken_x
bjorken_y     = options.do_bjorken_y
hadrons       = options.hadrons_cosz
energyQ       = options.energy
opening_angle = options.opening_angle
secondary_e   = options.secondary_e
hadrons_def   = options.hadrons_def
bjyve         = options.BjYvE

if no_weight:
    vscale = [10**-7, 10**0]
else:
    vscale = [10**3, 10**9]

print("Loading in the hdf5 files!")
import h5py

# open hdf5 file
op_file = h5py.File(ops_file,'r')
LepI_input = {}
# transcribe the contents into a dictionary. 
#   later on, we will access the contents of the keys and cast them as lists
for key in op_file:
    # sending the hdf5 list-like object to a ndarray and then to a list is WAY WAY WAY faster than directly casting it as a list
    LepI_input[key] = np.array(op_file[key]).tolist()
op_file.close()

ops_folder = "/".join(ops_file.split("/")[:-1])

# close file, open new file


# used to construct data_dict keys and make titles
GR_Mode = False
if not GR_Mode:
    flavors   = ['electron', 'muon', 'tau']
    matters   = ['matter', 'anti_matter']
    currents  = ['CC', 'NC'] #, 'GR']
else:
    flavors     = ['electron']
    matters     = ['anti_matter']
    currents    = ['GR']

class color_obj:
    def __init__(self, sep_flav, sep_mat, sep_cur):
        self.sep_flav = sep_flav
        self.sep_mat = sep_mat
        self.sep_cur = sep_cur
        
        self.n_colors = 1.0
        if self.sep_flav:
            self.n_colors *= len(flavors)
        if self.sep_mat:
            self.n_colors *= len(matters)
        if self.sep_cur:
            self.n_colors *= len(currents)

        self.cmap = plt.get_cmap( color )

    def get(self, flavor, matter, current):
        which_color = 0
        mult_scale = 1
        
        if self.sep_flav:
            which_color += mult_scale*flavor 
            mult_scale *= len(flavors)
        
        if self.sep_mat:
            which_color += mult_scale*matter
            mult_scale *= len(matters)

        if self.sep_cur:
            which_color += mult_scale*current

        colorVal = self.cmap( which_color/self.n_colors )
        assert( colorVal is not None )
        return( colorVal )

class title_obj:
    def __init__(self, sep_flav, sep_mat, sep_cur):
        assert( isinstance( sep_flav, bool))
        assert( isinstance( sep_mat, bool))
        assert( isinstance( sep_cur, bool ))

        self.sep_flav = sep_flav
        self.sep_mat  = sep_mat
        self.sep_cur  = sep_cur
        self.flavor = 0
        self.matter = 0
        self.current = 0

    def get(self, flavor, matter, current):
        assert( isinstance(flavor, int))
        assert( isinstance(matter,int))
        assert( isinstance( current, int))
        self.flavor = flavor
        self.matter = matter
        self.current = current

        title = ""
        if self.sep_flav:
            title += flavors[self.flavor]
        
        if self.sep_mat:
            if title!="":
                title+=", "
            
            if self.matter == 1:
                title += "anti-matter"
            else:
                title += matters[self.matter] 


        if self.sep_cur:
            if title!="":
                title+=", "
            title += currents[self.current] 

        return( title )

 
color = color_obj( sep_flavor, sep_matter, sep_current)        
title = title_obj( sep_flavor, sep_matter, sep_current) 

# function that constructs a title and file name for each plot
def make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current):
    title = ""
    if sep_flavor:
        if title=="":
            title += flavors[flavor]
        else:
            title += "_"+flavors[flavor]
    if sep_matter:
        if title=="":
            title += matters[matter]
        else:
            title += "_" + matters[matter]
    if sep_current:
        if title=="":
            title += currents[current]
        else:
            title += "_"+currents[current]
    return(title)


# will store data for plotting in a 3D array of lists
# first, need the size of each of these dimensions
flavD = 1
currD = 1
mattD = 1
if sep_flavor:
    flavD = len(flavors)
if sep_current:
    currD = len(currents)
if sep_matter:
    mattD = len(matters)
# ex - if all are false each dimension will only be 1 deep. 1x1x1 = 1 plot
#    - if just flavor is true, will have two 1 deep, 1 3 deep. 3x1x1 = 3 plots

# how many bins. This is used by like, all, of the plotting parts
nBins = 200

# should the Bjorken plots be log-scale? Nothing to do with writing a log
bjorken_log = True

# , x_err = None, y_err= None ,   

#get weights!
LepI_weights = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
use_weights= True
try:
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index. Still scanning over these dimension (to get the data), but putting it in the same plot data 3D list entry
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                
                LepI_weights[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_eventWeight"])
except KeyError:
    warnings.warn("Weight key not found, treating all as weight=1")
    use_weights = False

use_weights = use_weights and ( not no_weight )

plt.clf()

if bjorken_x:
    print("Starting process for Bjorken X") 
    
    # constuct the data storage objects from which to create the plots
    #plot data is stored in a 3D list, indexed by the "flavor index", "current index", and "matter index"
    #       earlier, we calculated the depth of these dimensions and stored them as 'currD' and the other '*D' things
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 


    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index. Still scanning over these dimension (to get the data), but putting it in the same plot data 3D list entry
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenX"])
                

    print("Making BJX Plots")
    # since the 3D list is now full of the data, we scan over each one and make its plot
    plt.figure()
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                if bjorken_log:
                    bins = np.logspace(-4,0,nBins)
                else:
                    bins = np.linspace( 0,1,nBins)

                # make as in calculated
                # make the LI histogram
                if use_weights:
                    histo_details  = np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), weights=LepI_weights[flavor][matter][current], bins=bins , density=False)[0]
                    scale_data  = histo_details/( bins[1:] - bins[:-1] )
                else:
                    histo_details = np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins )
                    scale_data = histo_details[0]/sum(histo_details[0])


                #plot the LI+NG histograms
                plt.errorbar(0.5*(bins[:-1]+bins[1:]) ,scale_data, color=color.get(flavor,matter,current), drawstyle='steps', label=title.get(flavor, matter, current) )

    # scale the axes, add labels, save fig
    plt.yscale('log')
    if bjorken_log:
        plt.xscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([vscale[0], vscale[1]])
    plt.legend()
    plt.ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
    plt.xlabel("Bjorken X", size=14)
    plt.savefig( ops_folder+"/output/BjorkenX.png", dpi=400)
    plt.clf()
else:
    print("Skipping Bjorken X") #do nothing...
plt.close()

# the same process is repeated for each one, with some slight changes for the words and axis ranges
# most of the remaining comments are similar to above, with exception forthe BjYvE plots
if bjorken_y:
    print("")
    print("")
    print("====> Starting process for Bjorken Y") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    tick_counter = 0
    print("Making BJY Plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                if bjorken_log:
                    bins = np.logspace(-4,0,nBins)
                else:
                    bins = np.linspace( 0,1,nBins)

                if use_weights:
                    histo_details  =  np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()),weights=LepI_weights[flavor][matter][current], bins=bins, density=False)[0]
                    scale_data = histo_details/( bins[1:]-bins[:-1])
                else:
                    # make the LI histogram
                    histo_details = np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])


                tick_counter += 1
                #plot the LI+NG histograms
                ax.errorbar(0.5*(bins[:-1]+bins[1:]) , scale_data, color=color.get(flavor,matter,current),ls='-', drawstyle='steps', label=title.get(flavor,matter,current))
                
    
    ax.set_yscale('log')
    if bjorken_log:
        ax.set_xscale('log')
    ax.set_xlim([bins.min(),bins.max()])
    ax.set_ylim([vscale[0],vscale[1]])
    ax.set_ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
    ax.legend()
    ax.set_xlabel("Bjorken Y", size=14)
    fig.savefig( ops_folder+"/output/"+ "BjorkenY_12.png", dpi=400)
    fig.clf()
    
    print("Tick counter at {}".format(tick_counter))
else:
    print("Skipping Bjorken Y")

plt.close()

if hadrons:
    print("")
    print("")
    print("====> Starting process for hadrons cos zenith angle") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_hadronsCosz"])

    print("")
    print("")
    print("====> Making Hadrons cos zenith plots")
    plt.figure()
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.linspace(-1,1,nBins)
                if use_weights:
                    histo_details  = np.histogram(np.cos(LepI_plots[flavor][matter][current]), weights=LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins, density=False)[0]
                    scale_data = histo_details/( bins[1:]-bins[:-1])
                else:
                    # make the LI histogram
                    histo_details = np.histogram(np.cos(LepI_plots[flavor][matter][current]),range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])

                #plot the LI+NG histograms
                plt.errorbar( 0.5*(bins[:-1]+bins[1:]) , scale_data ,color=color.get(flavor,matter,current),ls='-', drawstyle='steps', label=title.get(flavor,matter,current))
                
                #calculate the %difference and plot it

    plt.yscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([vscale[0], vscale[1]])
    plt.ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
    plt.legend()
    plt.xlabel(r"Hadrons $\cos\theta_{zenith}$", size=14)
    plt.savefig( ops_folder+"/output/"+ "HadronsCosZ.png", dpi=400)
    plt.clf()
else:
    print("Skipping Hadrons Cosz")
plt.close()

if opening_angle:
    print("")
    print("")
    print("====> Starting process for hadrons cos zenith angle") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_openingAngle"])
    plt.figure()
    print("")
    print("")
    print("====> Making opening angle plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.linspace(0,50,nBins)

                if use_weights:
                    histo_details  = np.histogram(np.array(LepI_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), weights=LepI_weights[flavor][matter][current], bins=bins, density=False)[0]
                    scale_data = histo_details/(bins[1:]-bins[:-1])
                else:
                    # make the LI histogram
                    histo_details = np.histogram(np.array(LepI_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])
                
                #plot the LI+NG histograms
                plt.errorbar( 0.5*(bins[:-1]+bins[1:]) , scale_data ,color=color.get(flavor,matter,current),ls='-', drawstyle='steps', label=title.get(flavor,matter,current))
                

    plt.yscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([vscale[0],vscale[1]])
    plt.legend()
    plt.ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
    plt.xlabel(r"Opening Angle: $\theta_{\mu}+\theta_{w}$", size=14)
    plt.savefig( ops_folder+"/output/"+ "openingAngle.png", dpi=400)
    plt.clf()
else:
    print("Skipping opening angle")
plt.close()

if secondary_e:
    print("")
    print("")
    print("====> Starting process for secondary lepton energy") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_secondaryE"])
    
    plt.figure()
    print("")
    print("")
    print("====> Making lepton energy plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.logspace(2,6,nBins)
                # make the LI histogram
                if use_weights:
                    histo_details  =  np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), weights=LepI_weights[flavor][matter][current], bins=bins, density=False)[0]
                    scale_data = 0.25*(bins[1:]+bins[:-1])*(bins[1:]+bins[:-1])*histo_details/(bins[1:]-bins[:-1])
                else:
                    histo_details = np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])


                #plot the LI+NG histograms
                plt.errorbar( 0.5*(bins[:-1]+bins[1:]) , scale_data ,color=color.get(flavor,matter,current),ls='-', drawstyle='steps', label=title.get(flavor,matter,current))
                
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([vscale[0], vscale[1]])
    plt.ylabel(r"$E^{2} dN/dE$ [GeV yr$^{-1}$]",size=14)
    plt.legend()
    plt.xlabel("Secondary Lepton Energy [GeV]", size=14)
    plt.savefig( ops_folder+"/output/"+ "secondaryE.png", dpi=400)
    plt.clf()
else:
    print("Skipping secondary energy")

plt.close()

if hadrons_def:
    print("")
    print("")
    print("====> Starting process for hadron deflection angle") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_hadronsDef"])
    plt.figure()
    print("")
    print("")
    print("====> Making hadron deflection plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.linspace(0,50,nBins)
                # make the LI histogram
                if use_weights:
                    histo_details  = np.histogram(np.array(LepI_plots[flavor][matter][current])*180./pi,weights=LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins, density=False)[0]
                    scale_data = histo_details/(bins[1:]-bins[:-1])
                else:
                    histo_details = np.histogram(np.array(LepI_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])

                #plot the LI+NG histograms
                plt.errorbar( 0.5*(bins[:-1]+bins[1:]) , scale_data ,color=color.get(flavor,matter,current),ls='-', drawstyle='steps', label=title.get(flavor,matter,current))
                
    plt.yscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([vscale[0],vscale[1]])
    plt.ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
    plt.legend()
    plt.xlabel("Hadrons deflection [deg]", size=14)
    plt.savefig( ops_folder+"/output/"+ "hadronDef.png", dpi=400)
    plt.clf()
else:
    print("Skipping hadron deflection")

plt.close()

if bjyve:
    print("")
    print("")
    print("====> Starting process for PrimaryE vs Bjorken Y plots") 
       
    # there is a seaborn plotter thing that does bascially what I'm doing, but it's suuuuuper slow
    nBinsss = 25

    # this function takes a list of pairs [ energy, Bjorken Y] and a list of bin edges
    def get_means(data, bins):
        holder = [[] for i in range(nBinsss) ] # create a list for each bin
        for entry in data: #scan over the pairs 
            bin_no = 0
            while entry[0] > bins[bin_no]: # find which energy bin the pair is associated with
                bin_no += 1
                if (bin_no)==len(bins):
                    break
            if not bin_no==len(bins):
                holder[bin_no].append(entry[1]) # and add the Bjorken Y from the pair to that bin's list
        
        means = []
        errors = []
        for entry in holder:
            means.append( np.mean(entry) )
            errors.append(( np.nanstd(entry))/np.sqrt(len(entry)) ) #std squared is variance

        return([means, errors]) # calculate the average of each bin-list's entries, return the mean and variance 

    if True:
        # constuct the data storage objects from which to create the plots
        LepI_energy = [[[[] for current in range(len(currents))] for matter in range(len(matters))] for flavor in range(len(flavors))] 
        LepI_bjy    = [[[[] for current in range(len(currents))] for matter in range(len(matters))] for flavor in range(len(flavors))] 

        # fill the data storage objects...
        print("Preparing Data")
        for current in range(len(currents)):
            if True:
                current_index = current
            else:
                current_index = 0
            for flavor in range(len(flavors)):
                if True:
                    flavor_index = flavor
                else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
                    flavor_index = 0
                for matter in range(len(matters)):
                    if True:
                        matter_index = matter
                    else:
                        matter_index = 0
                    # adding to the data storage things from which to make the plots
                    # scanning over all the separate keys in the data from pickles. 
                    
                    LepI_energy[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_primaryE"])
                    LepI_bjy[flavor_index][matter_index][current_index]    += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])
        


        flux_means = np.array([ np.array([0. for i in range(nBinsss)]), np.array([0. for i in range(nBinsss)])])

        for current in range(len(currents)):
            for flavor in range(len(flavors)):
                for matter in range(mattD):
                    LepI_plots =  [ LepI_energy[flavor][matter][current], LepI_bjy[flavor][matter][current] ]
                    LepI_plots = np.array(LepI_plots)
                    LepI_plots = np.transpose(LepI_plots)
                    bins = np.logspace( 2, 6, nBinsss)

                    LepI_means = get_means(LepI_plots, bins) 
                    # array( array(means), array(standard deviations) )
                    LepI_means = np.array([ np.array(LepI_means[0]), np.array(LepI_means[1]) ])
                    
                    if current==0 and flavor==1:
                        weight=15./40
                    else:
                        weight=1./40

                    flux_means[0] += (weight*LepI_means[0])
                    flux_means[1] += (weight*LepI_means[1])
                    

                    if flux_means[0][0] is np.nan:
                        print("{}, {}, {}".format(current, flavor, matter))
                        raise Exception("")

    plt.fill_between(bins, flux_means[0]-flux_means[1], flux_means[0]+flux_means[1], color='r', alpha=0.2, label="LI: Flux-averaged")

    # constuct the data storage objects from which to create the plots
    LepI_energy = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    LepI_bjy    = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    # fill the data storage objects...
    print("Preparing Data")
    for current in range(len(currents)):
        if sep_current:
            current_index = current
        else:
            current_index = 0
        for flavor in range(len(flavors)):
            if sep_flavor:
                flavor_index = flavor
            else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
                flavor_index = 0
            for matter in range(len(matters)):
                if sep_matter:
                    matter_index = matter
                else:
                    matter_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                
                LepI_energy[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_primaryE"])
                LepI_bjy[flavor_index][matter_index][current_index]    += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])
    
    print("")
    print("")
    print("====> Making BjY vs E plots")

    def get_nu(mattD):
        if mattD==0:
            return( r"$\nu$" )
        else:
            return( r"$\bar{\nu}$" )
    def get_curr(cur):
        if cur==0:
            return("CC")
        else:
            return("NC")
    

    print(mattD)
    plt.figure(1)
    


    for current in range(currD):
        #plt.clf()
        #figs, axes = plt.subplots(nrows=1, ncols=1, sharex=True, gridspec_kw={'height_ratios':[1]})
        for flavor in range(flavD): 
            for matter in range(mattD):
                # prepare the data
                if use_weights:
                    warnings.warn("Warning! I haven't figured out how to weight these ones properly...")
                print("Current: {}, Flavor: {}, Matter: {}".format(current, flavor, matter))
                try:
                    LepI_plots =  [ LepI_energy[flavor][matter][current], LepI_bjy[flavor][matter][current] ]
                    LepI_plots = np.array(LepI_plots)
                    LepI_plots = np.transpose(LepI_plots)
                except TypeError:
                    print( flavor )
                    print( matter )
                    print( current )
                    print( LepI_energy[flavor][matter][current] )
                    print( LepI_bjy[flavor][matter][current] )
                    raise Exception()
                
                # set up some bins, and then grab the mean Bjorken Y's we want to plot
                bins = np.logspace( 2, 6, nBinsss)

                LepI_means = get_means(LepI_plots, bins) 
                # array( array(means), array(standard deviations) )
                LepI_means = np.array([ LepI_means[0], LepI_means[1] ])
                
                # index zero means the data
                # index one means the error bars
                if current==0:
                    ls = '-'
                else:
                    ls ='--'
                #plt.errorbar(x=bins,y=LepI_means[0],yerr=LepI_means[1],xerr=None,capsize=5, drawstyle='steps', color=get_color(matter+2*current), label="{} ".format(currents[current])+get_nu(matter) )
                plt.errorbar(x=bins,y=LepI_means[0],yerr=None,xerr=None,capsize=5, drawstyle='steps', color=color.get(flavor,matter,current), label="LI: "+title.get(flavor,matter,current) )
#                axes[0].errorbar(x=bins,y=NuGe_means[0],yerr=NuGe_means[1],xerr=None,capsize=5, drawstyle='steps', color=get_color(2+matter), label="NG "+get_nu(matter) )
    


    # plot the data! 
    data_path = "/home/benito/Dropbox/Documents/Work/Research/IceCube/Projects/LeptonInjector/paper_data.json"
    data_obj = open( data_path, 'r')
    data = json.load( data_obj )
    data_obj.close()

    data_conv = []
    # <bad>
    for point in data:
        data_conv.append( np.array([data[point]["x"], data[point]["x_var"][0], data[point]["x_var"][1], data[point]["y"], data[point]["y_var"][0], data[point]["y_var"][1]]) )
    data_conv = np.array(data_conv).transpose()

    plt.hlines( data_conv[3], data_conv[1], data_conv[2], label="IceCube Inelasticity 5yr")
    plt.vlines( data_conv[0], data_conv[4], data_conv[5])


    # </bad>
       

    plt.xscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([0,0.6])
    plt.xlabel(r'$E_{\nu}$ [GeV]', size=16)
    plt.ylabel(r'$\left< y \right>$', size=16)            
    plt.legend(loc='lower left')
#    plt.grid(which='major', alpha=0.6)
#    plt.grid(which='minor', alpha=0.3)
    plt.tight_layout()
    plt.savefig( ops_folder+"/output/super_bjve.png", dpi=400 )
    plt.clf()
        #figs.clf()
        #plt.close(figs)
else:
    print("Skipping BjY vs E")


if energyQ:
    print("")
    print("")
    print("====> Starting process for energy plots") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 

    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_plots[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_primaryE"])

    plt.figure()
    print("")
    print("")
    print("====> Making primary energy plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.logspace(2,6,nBins)
                # make the LI histogram
                if use_weights:
                    histo_details  =  np.histogram(LepI_plots[flavor][matter][current],weights = LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins , density=False)[0]
                    scale_data = (((bins[:-1]+bins[1:])/2.)**2)*histo_details/(bins[1:]-bins[:-1])
                else:
                    histo_details = np.histogram(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins )
                    scale_data = histo_details[0]/sum(histo_details[0])
                
                plt.errorbar(0.5*(bins[:-1]+bins[1:]) ,scale_data,yerr=None,xerr=None,capsize=5, drawstyle='steps', color=color.get(flavor,matter,current), label="LI: "+title.get(flavor,matter,current) )
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([bins.min(),bins.max()])
    plt.ylim([vscale[0], vscale[1]])
    if no_weight:
        plt.ylabel(r"$E^{2}dN/dE$ [GeV]",size=14)
    else:
        plt.ylabel(r"$E^{2}dN/dE$ [GeV yr$^{-1}$]",size=14)
    plt.legend()
    plt.tight_layout()
    plt.xlabel("Primary Lepton Energy [GeV]", size=14)
    plt.savefig( ops_folder+"/output/"+ "primaryE.png", dpi=400)
    plt.clf()
else:
    print("Skipping primary energy")

if False:
    print("")
    print("")
    print("====> Starting process for Bjorken XvY Plots") 
    
    # constuct the data storage objects from which to create the plots
    LepI_bjx = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    LepI_bjy = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 


    print("Preparing Data")
    for flavor in range(len(flavors)):
        if sep_flavor:
            flavor_index = flavor
        else: # if we're not separating by flavor, the flavor-plot dimension is only 1 deep. Only use 0 index
            flavor_index = 0
        for matter in range(len(matters)):
            if sep_matter:
                matter_index = matter
            else:
                matter_index = 0
            for current in range(len(currents)):
                if sep_current:
                    current_index = current
                else:
                    current_index = 0
                # adding to the data storage things from which to make the plots
                # scanning over all the separate keys in the data from pickles. 
                LepI_bjx[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenX"])
                LepI_bjy[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])

    print("")
    print("")
    plt.close()
    print("====> Making Bj XvY Plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                # define the bins. These will be used for both the x and y axes
                nBins = 30
#                bins= np.linspace(10**-4, 1.0, nBins)
                bins = np.logspace(-4,0,nBins)
    
                # this is needed to put the colorbar on a log scale
                from matplotlib.colors import LogNorm
                # creat the plot title
              
                if use_weights:
                    normed_hist  = np.histogram2d(  LepI_bjx[flavor][matter][current] , LepI_bjy[flavor][matter][current], bins, weights=LepI_weights[flavor][matter][current])[0]
#                    normed_hist */ 
                else:
                    # calculate the Lepton Injector histogram
                    hist_details = np.histogram2d(  LepI_bjx[flavor][matter][current] , LepI_bjy[flavor][matter][current], bins )
                    # normalize the bins, and assign them a new name
                    normed_hist = np.array(hist_details[0])/len(LepI_bjx[flavor][matter][current])
                
                
                # use the matplotlib pcolormesh function to creat the histogram
                #      note - pcolormesh expects the data as (Y,X), and so I take the transpose of hte normalized histogram 
                
                # plt.pcolormesh( bins, bins, np.transpose(normed_hist), cmap = 'viridis', vmin=10**-6, vmax=1, norm=LogNorm())
                # cmap=color 
        
                

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                pos = 0.5*(bins[1:] + bins[:-1])
                wid = bins[1:]-bins[:-1]               
                x_data, y_data = np.meshgrid( pos, pos )
                widx, widy = np.meshgrid( wid, wid)
                
                normed_hist /= widx 
                normed_hist /= widy

                #x_data = x_data.flatten()
                #y_data = y_data.flatten()
                #widx = widx.flatten()
                #widy = widy.flatten()
                #z_data = normed_hist.flatten()
                
                #ax.bar3d( x_data, y_data, np.zeros(len(z_data)),1,1, z_data)
#                ax.bar3d( x_data, y_data, np.zeros(len(z_data)), widx, widy, z_data)


                # just add in all the dumb shit 
                ax.set_xlabel("log(Bjorken X)", size=18)
                ax.set_ylabel("log(Bjorken Y)", size=18)
                #ax.colorbar()
                #set(plt.gca(), 'colorscale', 'log')
                ax.set_xlim([np.log10(bins.min()), np.log10(bins.max())])
                ax.set_ylim([np.log10(bins.min()), np.log10(bins.max())])
#                ax.set_zscale('log')
#                ax.set_xscale('log')
#                ax.set_yscale('log')
                ax.set_zlabel(r"log(Flux $\frac{\partial N}{\partial X\partial Y}$) [yr$^{-1}$]", size=18)
                
                #ax.contourf3D( x_data,y_data, normed_hist, cmap ='viridis')
                suff = ax.plot_surface( np.log10(x_data), np.log10(y_data), np.log10(normed_hist), cmap='viridis',vmin=1,vmax=5, shade=True)
                #ax.title("Lepton Injector: " + title.get(flavor,matter,current), size=16)
                fig.colorbar(suff, shrink=0.5, aspect=5)
                plt.show()
#                plt.savefig( ops_folder+"/output/"+ "LeptonInjector_XvY"+title.get(flavor,matter,current)+".png", dpi=400)

else:
    print("Skipping Bjorken X vs Y plots")


