#!/cvmfs/icecube.opensciencegrid.org/py3-v4/RHEL_7_x86_64/bin/python3.6


# Ben Smithers
# benjamin.smithers@mavs.uta.edu



import h5py

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

# this is required for running as a condor job
os.environ['HDF5_USE_FILE_LOCKING'] ='FALSE'

pi = np.pi

parser = OptionParser()
parser.add_option("-i", "--ops", 
                  dest="ops_folder", 
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

force_unweighted = False

options, args = parser.parse_args()

# using lots of args to wrest control over how the plots are made
#store parser args! 
sep_flavor    = options.sep_flavor
sep_current   = options.sep_current
sep_matter    = options.sep_matter
ops_folder    = options.ops_folder
color         = options.color
bjorken_x     = options.do_bjorken_x
bjorken_y     = options.do_bjorken_y
hadrons       = options.hadrons_cosz
energyQ       = options.energy
opening_angle = options.opening_angle
secondary_e   = options.secondary_e
hadrons_def   = options.hadrons_def
bjyve         = options.BjYvE

lims = [10**-8, 10**0]

try:
    os.mkdir( ops_folder + "output" )
except FileExistsError: 
    pass

# prepare the color palette
this_cmap = plt.get_cmap(color)
def get_color(n, colormax=3.0):
    #0,1,2,3,
    colorVal = this_cmap(n/colormax)
    return(colorVal) # (r, g, b)



print("Loading in the hdf5 files!")
import h5py

# open hdf5 file
op_file = h5py.File(ops_folder + "LepI_processed.hdf5",'r')
LepI_input = {}
# transcribe the contents into a dictionary. 
#   later on, we will access the contents of the keys and cast them as lists
for key in op_file:
    # sending the hdf5 list-like object to a ndarray and then to a list is WAY WAY WAY faster than directly casting it as a list
    LepI_input[key] = np.array(op_file[key]).tolist()
op_file.close()

op_file = h5py.File(ops_folder + "NuGe_processed.hdf5",'r')
NuGe_input = {}
for key in op_file:
    NuGe_input[key] = np.array(op_file[key]).tolist()
op_file.close()


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
nBins = 400

# should the Bjorken plots be log-scale? Nothing to do with writing a log
bjorken_log = True

#get weights!
LepI_weights = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
NuGe_weights = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]
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
                NuGe_weights[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_eventWeight"])
except KeyError:
    warnings.warn("Weight key not found, treating all as weight=1")
    use_weights = False

use_weights = use_weights and (not force_unweighted)
if use_weights:
    the_lims = [10**-6, 10**-1]
else:
    the_lims = [10**-8, 10**0]

if bjorken_x:
    print("Starting process for Bjorken X") 
    
    # constuct the data storage objects from which to create the plots
    #plot data is stored in a 3D list, indexed by the "flavor index", "current index", and "matter index"
    #       earlier, we calculated the depth of these dimensions and stored them as 'currD' and the other '*D' things
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]


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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenX"])
                

    print("Making BJX Plots")
    # since the 3D list is now full of the data, we scan over each one and make its plot
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
                    histo_details  = plt.hist(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), weights=LepI_weights[flavor][matter][current], bins=bins, log=bjorken_log, density=False)[0]
                    histo_details2 = plt.hist(NuGe_plots[flavor][matter][current],range=(bins.min(),bins.max()), weights=NuGe_weights[flavor][matter][current], bins=bins, log=bjorken_log, density=False)[0]
                    scale_data  = histo_details/( bins[1:] - bins[:-1] )
                    scale_data2 = histo_details2/(bins[1:] - bins[:-1] )
                else:
                    histo_details = plt.hist(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins, log=bjorken_log)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(NuGe_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins,log=bjorken_log)
                    scale_data2 = histo_details[0]/sum(histo_details[0])
                
                # prepare the axes 
                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [(scale_data[i]-scale_data2[i])/scale_data[i] for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                # scale the axes, add labels, save fig
                axes[0].set_yscale('log')
                if bjorken_log:
                    axes[0].set_xscale('log')
                    axes[1].set_xscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                axes[1].set_ylim([-2,2])
                axes[0].legend()
                axes[0].set_ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel("Bjorken X", size=14)
                figs.savefig( ops_folder+"/output/"+ "BjorkenX_"+plot_title+".png", dpi=400)
                plt.close()
else:
    print("Skipping Bjorken X") #do nothing...


# the same process is repeated for each one, with some slight changes for the words and axis ranges
# most of the remaining comments are similar to above, with exception forthe BjYvE plots
if bjorken_y:
    print("")
    print("")
    print("====> Starting process for Bjorken Y") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])

    print("Making BJY Plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                if bjorken_log:
                    bins = np.logspace(-4,0,nBins)
                else:
                    bins = np.linspace( 0,1,nBins)

                if use_weights:
                    histo_details  =  plt.hist(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()),weights=LepI_weights[flavor][matter][current], bins=bins, log=bjorken_log, density=False)[0]
                    histo_details2 =  plt.hist(NuGe_plots[flavor][matter][current],range=(bins.min(),bins.max()),weights=NuGe_weights[flavor][matter][current], bins=bins, log=bjorken_log, density=False)[0]
                    scale_data = histo_details/( bins[1:]-bins[:-1])
                    scale_data2= histo_details2/(bins[1:]-bins[:-1])
                else:
                    # make the LI histogram
                    histo_details = plt.hist(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins, log=bjorken_log)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(NuGe_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins, log=bjorken_log)
                    scale_data2 = histo_details[0]/sum(histo_details[0])

                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [ (scale_data[i]-scale_data2[i])/scale_data[i]for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                axes[0].set_yscale('log')
                if bjorken_log:
                    axes[0].set_xscale('log')
                    axes[1].set_xscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                axes[1].set_ylim([-2,2])
                axes[0].set_ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
                axes[0].legend()
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel("Bjorken Y", size=14)
                figs.savefig( ops_folder+"/output/"+ "BjorkenY_"+plot_title+".png", dpi=400)
                plt.close()
else:
    print("Skipping Bjorken Y")

if hadrons:
    print("")
    print("")
    print("====> Starting process for hadrons cos zenith angle") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_hadronsCosz"])

    print("")
    print("")
    print("====> Making Hadrons cos zenith plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.linspace(-1,1,nBins)
                if use_weights:
                    histo_details  = plt.hist(np.cos(LepI_plots[flavor][matter][current]), weights=LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins, density=False)[0]
                    histo_details2 = plt.hist(np.cos(NuGe_plots[flavor][matter][current]), weights=NuGe_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins, density=False)[0]
                    scale_data = histo_details/( bins[1:]-bins[:-1])
                    scale_data2= histo_details2/(bins[1:]-bins[:-1])
                else:
                    # make the LI histogram
                    histo_details = plt.hist(np.cos(LepI_plots[flavor][matter][current]),range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(np.cos(NuGe_plots[flavor][matter][current]),range=(bins.min(),bins.max()), bins=bins)
                    scale_data2 = histo_details[0]/sum(histo_details[0])

                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [(scale_data[i]-scale_data2[i])/scale_data[i] for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                axes[0].set_yscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                axes[1].set_ylim([-2,2])
                axes[0].set_ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
                axes[0].legend()
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel(r"Hadrons $\cos\theta_{zenith}$", size=14)
                figs.savefig( ops_folder+"/output/"+ "HadronsCosZ_"+plot_title+".png", dpi=400)
                plt.close()
else:
    print("Skipping Hadrons Cosz")


if opening_angle:
    print("")
    print("")
    print("====> Starting process for hadrons cos zenith angle") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_openingAngle"])

    print("")
    print("")
    print("====> Making opening angle plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.linspace(0,50,nBins)

                if use_weights:
                    histo_details  = plt.hist(np.array(LepI_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), weights=LepI_weights[flavor][matter][current], bins=bins, density=False)[0]
                    histo_details2 = plt.hist(np.array(NuGe_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), weights=NuGe_weights[flavor][matter][current], bins=bins, density=False)[0]
                    scale_data = histo_details/(bins[1:]-bins[:-1])
                    scale_data2= histo_details2/(bins[1:]-bins[:-1])
                else:
                    # make the LI histogram
                    histo_details = plt.hist(np.array(LepI_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(np.array(NuGe_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), bins=bins)
                    scale_data2 = histo_details[0]/sum(histo_details[0])

                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [ (scale_data[i]-scale_data2[i])/scale_data[i] for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                axes[0].set_yscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                axes[1].set_ylim([-2,2])
                axes[0].legend()
                axes[0].set_ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel(r"Opening Angle: $\theta_{\mu}+\theta_{w}$", size=14)
                figs.savefig( ops_folder+"/output/"+ "openingAngle_"+plot_title+".png", dpi=400)
                plt.close()
else:
    print("Skipping opening angle")

if secondary_e:
    print("")
    print("")
    print("====> Starting process for secondary lepton energy") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_secondaryE"])

    print("")
    print("")
    print("====> Making lepton energy plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.logspace(2,6,nBins)
                # make the LI histogram
                if use_weights:
                    histo_details  =  plt.hist(LepI_plots[flavor][matter][current],weights=LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins,log=True, density=False)[0]
                    histo_details2 =  plt.hist(NuGe_plots[flavor][matter][current],weights=NuGe_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins,log=True, density=False)[0]
                    scale_data = histo_details/(bins[1:]-bins[:-1])
                    scale_data2=histo_details2/(bins[1:]-bins[:-1])
                else:
                    histo_details = plt.hist(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins,log=True)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(NuGe_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins,log=True)
                    scale_data2 = histo_details[0]/sum(histo_details[0])

                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [ (scale_data[i]-scale_data2[i])/scale_data[i] for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                axes[0].set_yscale('log')
                axes[0].set_xscale('log')
                axes[1].set_xscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                axes[0].set_ylabel(r"Flux $dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
                axes[1].set_ylim([-2,2])
                axes[0].legend()
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel("Secondary Lepton Energy [GeV]", size=14)
                figs.savefig( ops_folder+"/output/"+ "secondaryE_"+plot_title+".png", dpi=400)
                plt.close()
else:
    print("Skipping secondary energy")



if hadrons_def:
    print("")
    print("")
    print("====> Starting process for hadron deflection angle") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_hadronsDef"])

    print("")
    print("")
    print("====> Making hadron deflection plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.linspace(0,50,nBins)
                # make the LI histogram
                if use_weights:
                    histo_details  = plt.hist(np.array(LepI_plots[flavor][matter][current])*180./pi,weights=LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins, density=False)[0]
                    histo_details2 = plt.hist(np.array(NuGe_plots[flavor][matter][current])*180./pi,weights=NuGe_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins, density=False)[0]
                    scale_data = histo_details/(bins[1:]-bins[:-1])
                    scale_data2=histo_details2/(bins[1:]-bins[:-1])
                else:
                    histo_details = plt.hist(np.array(LepI_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), bins=bins)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(np.array(NuGe_plots[flavor][matter][current])*180./pi,range=(bins.min(),bins.max()), bins=bins)
                    scale_data2 = histo_details[0]/sum(histo_details[0])

                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [ (scale_data[i]-scale_data2[i])/scale_data[i] for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                axes[0].set_yscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                axes[0].set_ylabel(r"$dN/dE$ [GeV$^{-1}$yr$^{-1}$]",size=14)
                axes[1].set_ylim([-2,2])
                axes[0].legend()
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel("Hadrons deflection [deg]", size=14)
                figs.savefig( ops_folder+"/output/"+ "hadronDef_"+plot_title+".png", dpi=400)
                plt.close()
else:
    print("Skipping hadron deflection")

if bjyve:
    print("")
    print("")
    print("====> Starting process for PrimaryE vs Bjorken Y plots") 
       
    # there is a seaborn plotter thing that does bascially what I'm doing, but it's suuuuuper slow
    nBins = 50

    # this function takes a list of pairs [ energy, Bjorken Y] and a list of bin edges
    def get_means(data, bins):
        holder = [[] for i in range(nBins) ] # create a list for each bin
        for entry in data: #scan over the pairs 
            bin_no = 0
            while entry[0] > bins[bin_no]: # find which energy bin the pair is associated with
                bin_no += 1
                if bin_no==len(bins):
                    break
            if not bin_no==len(bins):
                holder[bin_no].append(entry[1]) # and add the Bjorken Y from the pair to that bin's list
        
        means = []
        errors = []
        for entry in holder:
            means.append( np.mean(entry) )
            errors.append(( np.nanstd(entry))/np.sqrt(len(entry)) ) #std squared is variance 
        return([means, errors]) # calculate the average of each bin-list's entries, return the mean and variance 

    print(currD)
    print(mattD)
    print(flavD)
    
    # constuct the data storage objects from which to create the plots
    LepI_energy = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_energy = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]
    LepI_bjy    = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_bjy    = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_energy[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_primaryE"])
                LepI_bjy[flavor_index][matter_index][current_index]    += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])
                NuGe_bjy[flavor_index][matter_index][current_index]    += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])
    
    for x in range(currD):
        for y in range(flavD):
            for z in range(mattD):
                print(type( LepI_energy[y][z][x]))
    print("flavors: {}".format(len(LepI_energy)))
    print("matters: {}".format(len(LepI_energy[0])))
    print("currents: {}".format(len(LepI_energy[0][0])))
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

    for current in range(currD):
        #plt.clf()
        figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
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
                NuGe_plots =  np.array([ NuGe_energy[flavor][matter][current], NuGe_bjy[flavor][matter][current] ]).transpose()
                
                # set up some bins, and then grab the mean Bjorken Y's we want to plot
                bins = np.logspace( 3, 6.5, nBins)
                LepI_means = get_means(LepI_plots, bins) 
                LepI_means = np.array([ LepI_means[0], LepI_means[1] ])
                NuGe_means = get_means(NuGe_plots, bins)
                NuGe_means = np.array([ NuGe_means[0], NuGe_means[1] ])
                
                # index zero means the data
                # index one means the error bars
                # yerr= LepI_means[1] and NuGe_means[1]
                axes[0].errorbar(x=bins,y=LepI_means[0],yerr=None,xerr=None,capsize=5, drawstyle='steps', color=get_color(0+matter), label="LI "+get_nu(matter) )
                axes[0].errorbar(x=bins,y=NuGe_means[0],yerr=None,xerr=None,capsize=5, drawstyle='steps', color=get_color(2+matter), label="NG "+get_nu(matter) )
                axes[0].set_xscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim([0,0.6])
                axes[1].set_xlabel(r'$E_{\nu}$ [GeV]', size=14)
                axes[0].set_ylabel('Bjorken Y', size=14)
                
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[1].set_xscale('log')
                axes[1].set_ylim([0.9,1.1])
                axes[1].errorbar(x=bins, y=LepI_means[0]/NuGe_means[0],yerr=(LepI_means[0]/NuGe_means[0])*(LepI_means[1]/LepI_means[0] + NuGe_means[1]/NuGe_means[0]),xerr=None,capsize=5, drawstyle='steps', color=get_color( 0.5+2*float(matter)), label="LI/NG "+get_nu(matter))
                

                #plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                #plt.title(plot_title, size=16)
                #plt.savefig( ops_folder+"/output/"+ "primaryE_vs_Y_"+plot_title+".png", dpi=400)
                #plt.close() 
        axes[0].legend()
        axes[1].legend()
        figs.savefig( ops_folder+"/output/super_bjve_{}.png".format(get_curr(current)), dpi=400 )
        figs.clf()
        #plt.close(figs)
else:
    print("Skipping BjY vs E")

if energyQ:
    print("")
    print("")
    print("====> Starting process for energy plots") 
    
    # constuct the data storage objects from which to create the plots
    LepI_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_plots = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]

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
                NuGe_plots[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_primaryE"])

    print("")
    print("")
    print("====> Making primary energy plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                bins = np.logspace(2,6,nBins)
                # make the LI histogram
                if use_weights:
                    histo_details  =  plt.hist(LepI_plots[flavor][matter][current],weights = LepI_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins,log=True, density=False)[0]
                    histo_details2 =  plt.hist(NuGe_plots[flavor][matter][current],weights = NuGe_weights[flavor][matter][current], range=(bins.min(),bins.max()), bins=bins,log=True, density=False)[0]
                    scale_data =  0.5*((bins[:-1]+bins[1:])**1)*histo_details/(bins[1:]-bins[:-1])
                    scale_data2=  0.5*((bins[:-1]+bins[1:])**1)*histo_details2/(bins[1:]-bins[:-1])
                else:
                    histo_details = plt.hist(LepI_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins,log=True)
                    scale_data = histo_details[0]/sum(histo_details[0])

                    # make the NG histogram
                    histo_details = plt.hist(NuGe_plots[flavor][matter][current],range=(bins.min(),bins.max()), bins=bins,log=True)
                    scale_data2 = histo_details[0]/sum(histo_details[0])

                plt.clf()
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
                figs, axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
                figs.suptitle(plot_title, size=16)

                #plot the LI+NG histograms
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data , color=get_color(0),ls='-', drawstyle='steps', label="Lepton Injector")
                axes[0].plot( 0.5*(bins[:-1]+bins[1:]) , scale_data2, color=get_color(3),ls='-', drawstyle='steps', label="NuGen DetMode Circle")
                
                #calculate the %difference and plot it
                # this is LI - NG
                ratios = [scale_data2[i]/scale_data[i] for i in range(len(scale_data))]
                axes[1].plot( 0.5*(bins[:-1]+bins[1:]) , ratios, color=(1,0,0),ls='-', drawstyle='steps')

                axes[0].set_yscale('log')
                axes[0].set_xscale('log')
                axes[1].set_xscale('log')
                axes[0].set_xlim([bins.min(),bins.max()])
                axes[1].set_xlim([bins.min(),bins.max()])
                axes[0].set_ylim(the_lims)
                #axes[1].set_ylim([0,0.5])
                axes[0].set_ylabel(r"$E^{1}dN/dE$ [yr$^{-1}$]",size=14)
                axes[0].legend()
                axes[1].set_ylabel("LI-NuG/LI", size=14)
                axes[1].set_xlabel("Primary Lepton Energy [GeV]", size=14)
                figs.savefig( ops_folder+"/output/"+ "primaryE_"+plot_title+".png", dpi=400)
                plt.close()
    Lepi_plots = None
    NuGe_plots = None
else:
    print("Skipping primary energy")

if False:
    print("")
    print("")
    print("====> Starting process for Bjorken XvY Plots") 
    
    # constuct the data storage objects from which to create the plots
    LepI_bjx = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_bjx = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]
    LepI_bjy = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)] 
    NuGe_bjy = [[[[] for current in range(currD)] for matter in range(mattD)] for flavor in range(flavD)]


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
                NuGe_bjx[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenX"])
                LepI_bjy[flavor_index][matter_index][current_index] += list(LepI_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])
                NuGe_bjy[flavor_index][matter_index][current_index] += list(NuGe_input[flavors[flavor]+"_"+matters[matter]+"_"+currents[current]+"_BjorkenY"])

    print("")
    print("")
    print("====> Making Bj XvY Plots")
    for flavor in range(flavD): 
        for matter in range(mattD):
            for current in range(currD):
                # define the bins. These will be used for both the x and y axes
                nBins = 20
                bins = np.logspace(-4,0,nBins)
    
                # this is needed to put the colorbar on a log scale
                from matplotlib.colors import LogNorm
                # creat the plot title
                plot_title = make_title(sep_flavor, sep_matter, sep_current, flavor, matter, current)
              
                if use_weights:
                    normed_hist  = np.histogram2d(  LepI_bjx[flavor][matter][current] , LepI_bjy[flavor][matter][current], bins, normed=True, weights=LepI_weights[flavor][matter][current])[0]
                    normed_hist2 = np.histogram2d(  NuGe_bjx[flavor][matter][current] , NuGe_bjy[flavor][matter][current], bins, normed=True, weights=NuGe_weights[flavor][matter][current])[0]
                else:
                    # calculate the Lepton Injector histogram
                    hist_details = np.histogram2d(  LepI_bjx[flavor][matter][current] , LepI_bjy[flavor][matter][current], bins )
                    # normalize the bins, and assign them a new name
                    normed_hist = np.array(hist_details[0])/len(LepI_bjx[flavor][matter][current])
                    
                    # do it again, but this time for NuGen. 
                    hist_details2 = np.histogram2d(  NuGe_bjx[flavor][matter][current] , NuGe_bjy[flavor][matter][current], bins )
                    
                    
                    normed_hist2 = np.array(hist_details2[0])/len(NuGe_bjx[flavor][matter][current])

                
                plt.clf() # clear any figure created by calling np.histogram2d
                # use the matplotlib pcolormesh function to creat the histogram
                #      note - pcolormesh expects the data as (Y,X), and so I take the transpose of hte normalized histogram 
                plt.pcolormesh( bins, bins, np.transpose(normed_hist), cmap = color, vmin=10**-6, vmax=1, norm=LogNorm())
                # cmap=color 

                # just add in all the dumb shit 
                plt.xlabel("Bjorken X", size=14)
                plt.xscale('Log')
                plt.ylabel("Bjorken Y", size=14)
                plt.yscale('Log')
                plt.colorbar()
                #set(plt.gca(), 'colorscale', 'log')
                plt.xlim([bins.min(), bins.max()])
                plt.ylim([bins.min(), bins.max()])
                plt.title("Lepton Injector: " + plot_title, size=16)
                plt.savefig( ops_folder+"/output/"+ "LeptonInjector_XvY"+plot_title+".png", dpi=400)
                plt.clf()

                
                plt.clf()
                plt.pcolormesh( bins, bins, np.transpose(normed_hist2), cmap = color, vmin=10**-6, vmax=1, norm=LogNorm())

                plt.xlabel("Bjorken X", size=14)
                plt.xscale('Log')
                plt.ylabel("Bjorken Y", size=14)
                plt.yscale('Log')
                plt.colorbar()
                #set(plt.gca(), 'colorscale', 'log')
                plt.xlim([bins.min(), bins.max()])
                plt.ylim([bins.min(), bins.max()])
                plt.title("NuGen: " + plot_title, size=16)
                plt.savefig( ops_folder+"/output/"+ "NuGen_XvY"+plot_title+".png", dpi=400)
                plt.clf()

                # do it again, again, but now look at the ratios between the two.
                ratios = (normed_hist/normed_hist2)
                plt.clf()
                plt.pcolormesh(bins, bins, np.transpose(ratios), cmap='coolwarm', vmin=0., vmax=2.)
                plt.xscale('Log')
                plt.yscale('Log')
                plt.xlabel('Bjorken X',size=14)
                plt.ylabel('Bjorken Y',size=14)
                plt.colorbar()
                plt.xlim([bins.min(), bins.max()])
                plt.ylim([bins.min(), bins.max()])
                plt.title("LepI / NuGe: "+ plot_title, size=16)
                plt.savefig( ops_folder+"/output/LI-NG-XvY"+plot_title+".png", dpi=400)
                plt.clf()
                plt.close()

else:
    print("Skipping Bjorken X vs Y plots")


