import GolemFitPy as gf # central fit utility

import json # used to load in a defaults file 
import os 

import numpy as np
from warnings import warn

# bring in a utility function for parsing the phase points 
from utils import parse_point, check_configuration

'''
Ben Smithers
benjamin.smithers@mavs.uta.edu

This script loads in a json file to configurea GolemFit Fit process.
A Realization file must be provided! 

It doesn't do anything yet, it's under active development 
'''


configuration_file = "configuration.json"
# using os to open the file correctly, regardless of how this script is called
f = open(os.path.join(os.path.dirname(__file__), configuration_file), 'r')
config = json.load(f)
f.close()

# load in configuration data
run_options= config['run_options']
check_configuration(run_options) # make sure there aren't any glaring errors
#       before wasting the user's time 
fit_config = config['central_values']
steering_config = config['steering_options']
parameters = config['parameters']
flags_config = config['fitflags']

# contruct GF objects
datapaths = gf.Datapaths(run_options['datapath'])
npp      = gf.NewPhysicsParams()
fitparams = gf.FitParameters(gf.sampleTag.Sterile)
steering_params = gf.SteeringParams(gf.sampleTag.Sterile)
fitparams_flag = gf.FitParametersFlag(True)
priors = gf.Priors(gf.sampleTag.Sterile)

# configure datapaths
point = run_options['point']
prompt_file = run_options['fluxdir']+'/prompt_atmospheric_'+point+'.hdf5'
astro_file  = run_options['fluxdir']+'/astro_'+point+'.hdf5'
conv_file   = run_options['fluxdir']+'/atmospheric_'+point+'.hdf5'
datapaths.conventional_nusquids_atmospheric_file = conv_file
datapaths.prompt_nusquids_atmospheric_file          = prompt_file
datapaths.astro_nusquids_file                       = astro_file
datapaths.barr_resources_location           = fluxdir

# Add configuration to GF objects 
for param in fit_config.keys():
    if not hasattr(fitparams, param)
        warn("fitparams doesn't have a parameter {}, but I tried setting it.".format(param))
    setattr(fitparams, param, fit_config[param]) 


for param in steering_config.keys():
    if not hasattr(steering_params, param):
        warn("steering params doesn't have a parameter {}, but I tried setting it".format(param))
    setattr(steering_params, param, steering_config[param])

for param in flags_config.keys():
    if not hasattr(fitparams_flag, param):
        warn("Fit Parameters Flag object doesn't have parameter {}, but I tried setting it".format(param))
    setattr(fitparams_flag, param, flags_config[param])

# set priors 

# set seeds 

golemfit = gf.GolemFit(datapaths, steering_params, npp)
golemfit.SetFitParametersFlag(fitparams_flag)
golemfit.SetFitParametersPriors(priors)


realization_dist = np.load(run_options['realization'])['realization']


