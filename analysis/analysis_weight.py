#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/snobo/py3-v4.0.1/


##!/usr/bin/python
##!/cvmfs/icecube.opensciencegrid.org/py3-v4/RHEL_7_x86_64/bin/python

# this processes the I3 file from a generation using LeptonWeighter or NuGen
# it also calculates the weights of each event according to which generation method was used

# combines the functionality of the 'analysis.py' and the 'process-Weight.py' scripts

# === Benjamin Smithers 
# === benjamin.smithers@mavs.uta.edu



#import non-icecube dependencies

# used for parsing inputs
from optparse import OptionParser
# used for verifying files exist
import os
import numpy as np
import h5py
from glob import glob

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
pi = np.pi
proton_mass = 0.93827231 #GeV


from I3Tray import I3Tray
from icecube import icetray, dataio, dataclasses
from icecube.weighting import weighting, get_weighted_primary
from icecube.icetray import I3Units
from icecube import LeptonInjector
import LeptonWeighter as LW



#specify input parameters
parser = OptionParser()
parser.add_option("-t", "--NuGen",
                  dest='NuGen',
                  default=False,
                  action="store_true",
                  help="Signal this as a NuGen hdf5 file.")
parser.add_option("-i","--input",
                  dest="input",
                  default="/data/user/bsmithers/runs/snobo_danger/processed/",
                  help="Where to find the input files")
options, args = parser.parse_args()

GR = True
print("Weighting to: " + ("GR" if GR else "NC/CC"))

#out_folder   = "/data/user/bsmithers/runs/grad" 
#out_folder = "/data/user/bsmithers/runs/snobo_danger/processed/"
in_folder = options.input
out_folder = in_folder+"/processed/"
NuGen        = options.NuGen
calc_weights = True

if not os.path.exists(out_folder):
    os.mkdir(out_folder)

# nugen file
#in_files = ["/data/user/bsmithers/runs/grad/nugen_gen.i3.zst"]

# LI data
in_files = glob(in_folder+"/*.i3.zst")
lic_files = glob(in_folder+"/*.lic")

print("Running over {} files".format(len(in_files)))

# not using the poly ones 
def pop_matching(skip_vals, in_list): 
    index = 0

    while index<len(in_list):
        if skip_vals in in_list[index]:
            in_list.pop(index)
        else:
            index += 1
    return(in_list) 

skip_vals = "_poly"
in_files = pop_matching(skip_vals, in_files)
lic_files = pop_matching(skip_vals, lic_files)

flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
liveTime   =3.1536e7

if not calc_weights:
    def get_weight(frame):
        return(1.)
else:
    if NuGen:
        print("Using NuGen format")
        def get_weight(frame):
            nGenerated  = frame['I3MCWeightDict']['NEvents']
            oneweight   = frame['I3MCWeightDict']['OneWeight']
            weight      = oneweight * liveTime / nGenerated 
            flux        = flux_params['constant']*( frame['NuGPrimary'].energy / flux_params['scale'])**flux_params['index']
            return( weight*flux  )
    else:


        print("Using LeptonInjector format")
        
        net_generation          = []
        for lic in lic_files:
            net_generation += LW.MakeGeneratorsFromLICFile( lic )
        print("finished building the generators")
        
        # assuming physical flux is a powerlaw obeying E^-2
        flux                    = LW.PowerLawFlux(flux_params['constant'] , flux_params['index'] , flux_params['scale'] )    

        # cross section follows 
        if GR:
            xs = LW.GlashowResonanceCrossSection()
        else:
        # these are the same ones in the cvmfs, just copied to my directory 
            xs = LW.CrossSectionFromSpline(
                    "/data/user/bsmithers/cross_sections/dsdxdy_nu_CC_iso.fits",
                    "/data/user/bsmithers/cross_sections/dsdxdy_nubar_CC_iso.fits",
                    "/data/user/bsmithers/cross_sections/dsdxdy_nu_NC_iso.fits",
                    "/data/user/bsmithers/cross_sections/dsdxdy_nubar_NC_iso.fits")


        weight_event 			= LW.Weighter(flux, xs, net_generation)

        def get_weight(fame):
            LWevent = LW.Event()
            EventProperties 				= frame['EventProperties']
            LeptonInjectorProperties 		= frame['LeptonInjectorProperties']
            LWevent.primary_type 			= LW.ParticleType(EventProperties.initialType)
            LWevent.final_state_particle_0 	= LW.ParticleType(EventProperties.finalType1)
            LWevent.final_state_particle_1 	= LW.ParticleType(EventProperties.finalType2)
            LWevent.zenith 					= EventProperties.zenith
            LWevent.energy 					= EventProperties.totalEnergy
            LWevent.azimuth 				= EventProperties.azimuth
            LWevent.interaction_x 			= EventProperties.finalStateX
            LWevent.interaction_y 			= EventProperties.finalStateY
            if isinstance(EventProperties, LeptonInjector.VolumeEventProperties):
                LWevent.radius 					= EventProperties.radius
                LWevent.total_column_depth      = 0.
            else:
                LWevent.total_column_depth 		= EventProperties.totalColumnDepth
                LWevent.radius                  = EventProperties.impactParameter
            LWevent.x                       = 0. 
            LWevent.y                       = 0.
            LWevent.z                       = 0.
            weight                          = weight_event( LWevent ) #LW function defined above! 
            return( weight*liveTime )


print("Starting!")

# returns the angle between two unit vectors defined by angles in spherical coordinates
#   slow, but it works 
def angle_between( theta1, phi1, theta2, phi2):


    x1 = np.cos(phi1)*np.sin(theta1)
    y1 = np.sin(phi1)*np.sin(theta1)
    z1 = np.cos(theta1)

    x2 = np.cos(phi2)*np.sin(theta2)
    y2 = np.sin(phi2)*np.sin(theta2)
    z2 = np.cos(theta2)

    cosangle = x1*x2 + y1*y2 + z1*z2
    return(np.arccos( cosangle )) #arg returned in radians! 

#returns true if the particle is matter
#returns false if it's antimatter
#       this is the dumbest way I could've implemented this, but it works 
#       ... for leptons only
def is_matter(particle):
    pType = str(particle.type)

    # NuE is less than four letters long, and so the second `if` statement doesn't work.
    #       so if the type is fewer than four letters, it's a NuE and therefore matter
    if len(pType)<4:
        return("matter")
    # this is even dumber 


    # all the antimatter particles end with either `Bar' or `Plus'
    if (pType[-3:]=='Bar' or pType[-4]=='Plus'):
        return("anti_matter")
    else:
        return("matter")

# just, look at the type. Is it one of the six neutrinos?
def is_neutrino(particle):
    pType = str(particle.type)
    if (pType == 'NuE' or pType == 'NuEBar' 
                or pType == 'NuMu' or pType == 'NuMuBar'
                or pType == 'NuTau' or pType == 'NuTauBar'):
        return(True)
    else:
        return(False)

def is_charged_lepton(particle):
    #if the particle is a lepton, and not a neutrino, return True
    #   also returns true for anti-leptons
    pType = str(particle.type)
    if (pType == 'EPlus' or pType == 'EMinus' 
                or pType == 'MuPlus' or pType == 'MuMinus'
                or pType == 'TauPlus' or pType == 'TauMinus'):
        return(True)
    else:
        return(False)

def get_flavor( particle ):
    pType = str(particle.type)
    E_flavor  = ['EPlus', 'EMinus', 'NuE', 'NuEBar']
    Mu_flavor = ['MuPlus', 'MuMinus', 'NuMu', 'NuMuBar']
    Tau_flavor= ['TauPlus', 'TauMinus', 'NuTau', 'NuTauBar']
    
    if pType in E_flavor:
        return('electron')
    elif pType in Mu_flavor:
        return('muon')
    elif pType in Tau_flavor:
        return('tau')
    else:
        # this will obviously need to change in some bsm searches
        raise Exception("That's not a lepton?")

GR_Mode = False
if GR_Mode:
    print("Cutting all non-GR events")
    def keep_event(interaction, primary, secondary, hadrons):
        if interaction =='GR':
            return(True)
        else:
            return(False)
else:
    def keep_event(key_prefix, primary, secondary,hadrons):
        return(True)

# which keys will be used
if not GR_Mode:
    all_flavors     = ['electron', 'muon', 'tau']
    all_matter      = ['matter', 'anti_matter']
    all_interactions= ['GR', 'NC', 'CC']
else:
    all_flavors     = ['electron']
    all_matter      = ['anti_matter']
    all_interactions= ['GR']


from collections import deque

# deques have O(n) insertion of elements , bit better than lists
data_dict = { }
dtypes = ['primaryE', 'openingAngle',  'secondaryE', 'BjorkenX', 'BjorkenY', 'hadronsCosz', 'hadronsDef', 'eventWeight']
for flavor in all_flavors:
    for matt in all_matter:
        for interaction in all_interactions:
            for dtype in dtypes:
                data_dict[flavor+"_"+matt+"_"+interaction+"_"+dtype] = deque([])


# stores the data in the data_dict object
def store_data(key_prefix, primary, secondary, hadrons, event_weight): 
    data_dict[key_prefix+"_primaryE"].append( primary.energy) 
    data_dict[key_prefix+"_secondaryE"].append( secondary.energy )
    data_dict[key_prefix+"_openingAngle"].append(  angle_between( 
            hadrons.dir.zenith, hadrons.dir.azimuth, 
            secondary.dir.zenith, secondary.dir.azimuth ) 
            )
    data_dict[key_prefix+"_hadronsCosz"].append(hadrons.dir.zenith )
    data_dict[key_prefix+"_hadronsDef"].append( angle_between( 
            hadrons.dir.zenith, 
            hadrons.dir.azimuth, 
            primary.dir.zenith, 
            primary.dir.azimuth) 
            )
    data_dict[key_prefix+"_BjorkenY"].append( 1.0 - (secondary.energy/primary.energy) )
    nu            = primary.energy - secondary.energy
    scatter_theta =  angle_between( 
            secondary.dir.zenith, 
            secondary.dir.azimuth, 
            primary.dir.zenith, 
            primary.dir.azimuth )
    sin_theta2    = np.sin( scatter_theta/2.0)**2
    Q_squared     = 4.*primary.energy*secondary.energy*sin_theta2
    data_dict[key_prefix+"_BjorkenX"].append(Q_squared / (2.*proton_mass*nu) )
    if isinstance(event_weight, float) or isinstance(event_weight,int):
        data_dict[key_prefix+"_eventWeight"].append( event_weight)
    else:
        data_dict[key_prefix+"_eventWeight"].append( event_weight.value )


#
#
# V-- this is where the actual bulk of the code is --V
#
#

for input_file in in_files:

    if not os.path.isfile(input_file):
        print("Uh oh! This file doesn't exist: "+ input_file )
        print("Skipping.")
        continue
    else:
        print("Now loading... "+ input_file)

    data_file = dataio.I3File( input_file, 'r')

    # scan over the frames
    while data_file.more():
        frame = data_file.pop_frame()
        if str(frame.Stop) != 'DAQ':
            #if this isn't a data frame, skip!
            continue

        # get the particles
        particle_list = list(frame['I3MCTree'])
        
        if len(particle_list)!=3:
            # if for some reason there isn't the right number of particles, skip! 
            print("There are more than three particles. This is not good. ")
            continue
        
        primary = None  # neutrino which caused the event
        hadrons = None  # hadronic shower produced
        secondary = None  # secondary nu/l produced

        interaction = "" 
        
        # how do we get the interaction types?? 
        if NuGen:
            interactiont = int(frame['I3MCWeightDict']['InteractionType'])
            if interactiont == 1:
                interaction = 'CC'
            elif interactiont == 2:
                interaction = 'NC'
            elif interactiont == 3:
                interaction = 'GR'
            else: 
                print("error?")
                continue
        else:
            clepton_count   = 0
            neutroh_count   = 0
            hadron_count    = 0
            for particle in particle_list: 
                if str(particle.shape)!='Primary':
                    if is_neutrino(particle):
                        neutroh_count+=1 # if it's a neutrino
                    elif is_charged_lepton(particle):
                        clepton_count+=1 # not neut, but lepton --> charged lepton
                    else:
                        hadron_count+=1  # otherwise it's 'hadrons'
                
            if hadron_count==1:  # 1 'hadrons' means NC or CC
                if clepton_count==1:
                    interaction='CC'  # charged...
                else:
                    interaction='NC' # not charged
            else:
                interaction='GR'
        
        
        if interaction!='GR':
            # use the CC/NC particle assignment algorithm
            for particle in particle_list:
                if str(particle.shape)=='Primary':
                    primary     = particle
                elif (is_charged_lepton(particle)) or (is_neutrino(particle)):
                    secondary   = particle
                else:
                    hadrons     = particle
        else: 
            # use the GR algorithm
            for particle in particle_list:
                # if it's the primary, make it the primary
                if str(particle.shape)=='Primary':
                    primary = particle
                # otherwise, check if 'secondary' is assigned
                #       if not make this particle the 'secondary'
                #       if 'secondary' is assigned, then we're at the last one and so call it 'hadrons'
                #           even though it really might just be leptons. That's just a name. Deal with it. 
                else:
                    if secondary is None:
                        secondary = particle
                    else:
                        hadrons   = particle


        #ideally, all three would be assigned. 
        #   if one of them isn't, let the user know and move on to the next particle

        if ( (primary is None) or (secondary is None) or (hadrons is None)):
            print("Something Went Wrong!")
            continue
        
        
        key = get_flavor(primary) +"_" + is_matter(primary)+"_"+interaction
        # call the store data function with the desired key
        if keep_event(interaction , primary, secondary, hadrons):
            store_data(key, primary, secondary, hadrons, get_weight(frame) )

    data_file.close()

print("Saving File")
if NuGen:
    name = "NuGe"
else:
    name = "LepI"

file_name = "{}/{}_processed.hdf5".format( out_folder, name) 


# open the file. Save the keys+data inside
phil = h5py.File(file_name, 'w')
for key in data_dict:
    dset = phil.create_dataset(key, data = data_dict[key])
phil.close()
