#!/cvmfs/icecube.opensciencegrid.org/py3-v4/RHEL_7_x86_64/bin/python

# === Benjamin Smithers 
# === benjamin.smithers@mavs.uta.edu


"""
Does some analysis on and calculates weights for events from a standalone lepton injector run

See:
    input_files
    Flux
    cross sections
"""


#import non-icecube dependencies

# used for verifying files exist
import os
import numpy as np
import sys

import LeptonWeighter as LW
# import h5py as h5


out_folder = "/home/benito/Desktop/testing/with_volume/"
#input_files = ["/home/bsmithers/processing/LI_standalone_test/data_output_0a.h5"]
input_files = ["/home/benito/Desktop/testing/data_output_0a.h5"]
#inp = "/home/bsmithers/processing/LI_standalone_test/data_output_0a.h5"

lic_files = ["/home/benito/Desktop/testing/config_0a.lic"]
#lic_files = ["/home/bsmithers/processing/LI_standalone_test/config_0a.lic"]

name = "LepI"

file_name = "{}/{}_processed.hdf5".format( out_folder, name) 

# outputs to "LepI_processed.hdf5" normally
# if this is set to True, it always does that
# otherwise it will append a number
overwrite = False

#input_files = glob("./data_output*.h5")


pi = np.pi
proton_mass = 0.93827231 #GeV

#Flux = '/home/bsmithers/software_dev/Analysis/nusquids_stuff/'
Flux = '/home/carguelles/vault/golem_fit_installation_test/sources/GolemFit/resources/Fluxes/conventional/'# conventional_atmospheric.hdf5'
#Flux = '/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/Fluxes/Flux_AIRS_sib_HG_th24_dm2/prompt_atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5'
#Flux = '/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/Resources/FluxOscCalculator_US/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5'
#Flux = '/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/Resources/GolemFit/resources/FluxOscCalculator/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5'
#Flux = '/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/Resources/FluxOscCalculator_US/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5'

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
    if particle>0:
        return("matter")
    else:
        return("anti_matter")

def get_flavor( particle ):
    numb = abs(particle)
    if numb==11 or numb==12:
        return("electron")
    elif numb==13 or numb==14:
        return("muon")
    elif numb==15 or numb==16:
        return("tau")
    else:
        raise ValueError("NOT A LEPTON")

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

if not GR_Mode:
    all_flavors     = ['electron', 'muon', 'tau']
    all_matter      = ['matter', 'anti_matter']
    all_interactions= ['GR', 'NC', 'CC']
else:
    all_flavors     = ['electron']
    all_matter      = ['anti_matter']
    all_interactions= ['GR']


data_dict = { }
dtypes = ['primaryE', 'openingAngle',  'secondaryE', 'BjorkenX', 'BjorkenY', 'hadronsCosz', 'hadronsDef', 'eventWeight']
for flavor in all_flavors:
    for matt in all_matter:
        for interaction in all_interactions:
            for dtype in dtypes:
                data_dict[flavor+"_"+matt+"_"+interaction+"_"+dtype] = []

flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
liveTime   =3.1536e7


net_generation = []
for lic in lic_files:
    net_generation += LW.MakeGeneratorsFromLICFile(lic) 

xs_folder = "/home/benito/Desktop/testing/resources/"
xs = LW.CrossSectionFromSpline(
                    os.path.join(xs_folder,"dsdxdy_nu_CC_iso.fits"),
                    os.path.join(xs_folder,"dsdxdy_nubar_CC_iso.fits"),
                    os.path.join(xs_folder,"dsdxdy_nu_NC_iso.fits"),
                    os.path.join(xs_folder,"dsdxdy_nubar_NC_iso.fits"))

use_atm = False
if use_atm:
    # commenting these out while I make a BSM flux 
    pi_nusquids_flux = LW.nuSQUIDSAtmFlux(Flux+'/pion_atmospheric.hdf5')
    k_nusquids_flux = LW.nuSQUIDSAtmFlux(Flux +'/kaon_atmospheric.hdf5')
    weight_pi = LW.Weighter( pi_nusquids_flux, xs, net_generation )
    weight_k = LW.Weighter( k_nusquids_flux, xs, net_generation )
    
    #the_flux = LW.nuSQUIDSAtmFlux( Flux + "fluxes_flavor.hdf5" )
    #weight_it = LW.Weighter( the_flux, xs, net_generation )

else:
    flux = LW.PowerLawFlux( flux_params['constant'] , flux_params['index'] , flux_params['scale'] )

    weight_event = LW.Weighter( flux, xs, net_generation )

import h5py as h5



data_file = h5.File( inp, 'r')
which = list(data_file.keys())
print(data_file[which[0]])
print(data_file[which[0]]['properties'])
data_file.close()

print("Passed the weird test")

called_once = False
def get_weight( props ):
    """
    Accepts the properties list of an event and returns the weight
    """
    LWevent = LW.Event()
    LWevent.energy = props[0]
    LWevent.zenith = props[1]
    LWevent.azimuth = props[2]
    
    LWevent.interaction_x = props[3]
    LWevent.interaction_y = props[4]
    LWevent.final_state_particle_0 = LW.ParticleType( props[5] )
    LWevent.final_state_particle_1 = LW.ParticleType( props[6] )
    LWevent.primary_type = LW.ParticleType( props[7] )
    LWevent.radius = props[8]
    LWevent.total_column_depth = props[9]
    LWevent.x = 0
    LWevent.y = 0
    LWevent.z = 0
    global called_once 
    if use_atm:
        try:
            #one_weight = weight_it( LWevent )
            one_weight = (weight_pi(LWevent) + weight_k(LWevent))*0.5
            if not called_once:
                called_once = True
                print("Called at least once!")
        except Exception:
            one_weight = 0.0
    else:
        one_weight = weight_event(LWevent)

    if one_weight==np.nan:
        raise ValueError("Bad Weight!")
    return( one_weight*liveTime )

# stores the data in the data_dict object
def store_data(key_prefix, primary, secondary, hadrons, this): 
    data_dict[key_prefix+"_primaryE"].append( primary[-1]) 
    data_dict[key_prefix+"_secondaryE"].append( secondary[-1] )
    data_dict[key_prefix+"_openingAngle"].append(  angle_between( 
            hadrons[3][0], hadrons[3][1], 
            secondary[3][0], secondary[3][1] ) 
            )
    data_dict[key_prefix+"_hadronsCosz"].append(hadrons[3][0] )
    data_dict[key_prefix+"_hadronsDef"].append( angle_between( 
            hadrons[3][0], 
            hadrons[3][1], 
            primary[3][0], 
            primary[3][1]) 
            )
    data_dict[key_prefix+"_BjorkenY"].append( 1.0 - (secondary[-1]/primary[-1]) )
    nu            = primary[-1] - secondary[-1]
    scatter_theta =  angle_between( 
            secondary[3][0], 
            secondary[3][1], 
            primary[3][0], 
            primary[3][1] )
    sin_theta2    = np.sin( scatter_theta/2.0)**2
    Q_squared     = 4.*primary[-1]*secondary[-1]*sin_theta2
    data_dict[key_prefix+"_BjorkenX"].append(Q_squared / (2.*proton_mass*nu) )
    data_dict[key_prefix+"_eventWeight"].append(this)

#
#
# V-- this is where the actual bulk of the code is --V
#
#

for input_file in input_files: 

    if not os.path.isfile(input_file):
        print("Uh oh! This file doesn't exist: "+ input_file )
        print("Skipping.")
        raise Exception("SKIP")
    else:
        print("Now loading... "+ input_file)

    data_file = h5.File( inp , 'r')
    key_list = list( data_file.keys() )

    # scan over the frames
    for injector in key_list:
        print("     + "+injector)
        number = len(np.array(data_file[injector]['properties'] ))
        for event in range(number):

            primary = data_file[injector]['initial'][event]  # neutrino which caused the event
            hadrons = data_file[injector]['final_2'][event]  # hadronic shower produced
            secondary = data_file[injector]['final_1'][event]  # secondary nu/l produced
            props = data_file[injector]['properties'][event]

            interaction = '' 
            
            if hadrons[1] == -2000001006:
                if secondary[1] in [ 12, 14, 16, -12, -14, -16]: 
                    interaction = 'NC'
                elif secondary[1] in [ 11, 13, 15, -11, -13, -15]:
                    interaction = 'CC'
                else:
                    raise Exception("Unexpected pdg id {}".format( secondary [1] ))
            else:
                if (hadrons[1]==-12 and secondary[1]==11) or (hadrons[1]==-14 and secondary[1]==13) or (hadrons[1]==-16 and secondary[1]==15) or (hadrons[1]==secondary[1]==-2000001006):
                    interaction='GR'
                else:
                    print("Not sure what to do with these: {} and {}".format(hadrons[1], secondary[1]))
                    continue
            if interaction!='': 
                key = get_flavor(secondary[1]) +"_" + is_matter(secondary[1])+"_"+interaction
                # call the store data function with the desired key
                if keep_event(interaction , primary, secondary, hadrons):
                    this = get_weight(props)

                    store_data(key, primary, secondary, hadrons, this)
            else:
                print("Somehow the int is still ''")

    data_file.close()

print("Saving File")

# Write the file! 


# this guy checks if I already have something with the same name, adds a -1 if so
# if there is already a numbered file by the same name, it makes the number bigger! 

while (os.path.exists(file_name) and (not overwrite)):
    # we have a file name

    # there may be periods in folder names, 
    # so we get all the possible early-splittings and rejoin them
    prefix = file_name.split(".")[:-1]
    prefix = ".".join(prefix)

    suffix = file_name.split(".")[-1]

    # check, is that last char a number?
    # there's probably a less 'plan on failing' method, but this one was easy to implement! 
    try:
        prefix_number = int(prefix[-1])
        # increase if so increase the number
        prefix_number += 1
        
        # and replace that last character with the new number!
        prefix = prefix[:-1] + str( prefix_number )
        
        if prefix_number>(1e6):
            raise Exception("You have too many files: {}".format(prefix_number)) # wait why did I... oh never mind
        
    except ValueError:
        #it's not a number at the end, 
        #   so we stick a number at the end!
        prefix = prefix + "-1"
    
    # and we have a new file name which either now has a "-1" or a number one higher
    file_name = ".".join([prefix, suffix])


# open the file. Save the keys+data inside
phil = h5.File(file_name, 'w')
for key in data_dict:
    dset = phil.create_dataset(key, data = data_dict[key])
phil.close()
