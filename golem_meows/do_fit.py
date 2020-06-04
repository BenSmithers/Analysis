import GolemFitPy as gf # central fit utility

import json # used to load in a defaults file 
import os # used to load in the configuration file

import numpy as np

# bring in several utility funtions 
from utils import parse_point, check_configuration, converter, get_seed, set_GF, explicit_convert
from utils import listattr

'''
Ben Smithers
benjamin.smithers@mavs.uta.edu

The script does the fit. Hence the name. DO FIT. It DOES THE FIT. 

This script loads in a json file to configure a GolemFit Fit process.
A Realization file must be provided! 

This script passes a number of flags, priors, and various other parameters to a GolemFit object.
It then does a LLH minimization with GF.

Then it outputs a set of parameters from the fit to a json file.
'''


configuration_file = "configuration.json"
# using os to open the file correctly, regardless of how this script is called
f = open(os.path.join(os.path.dirname(__file__), configuration_file), 'r')
config = json.load(f)
f.close()

# load in configuration data
# do the implicit conversion to make sure all the strings are strings! 
run_options= explicit_convert(config['run_options'])
# verify that these run options are valid before doing any heavy lifting
check_configuration(run_options)

fit_config      = explicit_convert(config['central_values'])
steering_config = explicit_convert(config['steering_options'])
parameters      = explicit_convert(config['parameters'])
flags_config    = explicit_convert(config['fitflags'])
priors_config   = explicit_convert(config['priors'])

# contruct GF objects
datapaths = gf.DataPaths(run_options['datapath'])
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
#fastMC_location = '/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/FastMC/'+ sim+'_'+fast_scaling+'_'+years+'_'+fluxdir.split('/')[-2]+'/'+ point+'/'

compact_file_path = '/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/FastMC/' + \
        steering_config['simToLoad'] +'_' + "{:.2f}".format(steering_config['fastmode_scaling']) +'_' + \
        str(parameters['years']) + '_' +run_options['fluxdir'].split('/')[-2] + \
        '/' + run_options['point'] + '/'
datapaths.compact_file_path = compact_file_path 
datapaths.conventional_nusquids_atmospheric_file = conv_file
datapaths.prompt_nusquids_atmospheric_file       = prompt_file
datapaths.astro_nusquids_file                    = astro_file
datapaths.barr_resources_location               = run_options['fluxdir']

# This Loads the parameters from the json file and injects them into the GF objects 
print("Setting the parameters")
set_GF(fitparams, fit_config)
set_GF(steering_params, steering_config) #previously caused the segault
set_GF(fitparams_flag, flags_config)
set_GF(priors, priors_config)
set_GF(npp, parse_point(point))

# Set seeds 

# set special settings
def steer():
    """
    This sets up special settings for the steering parameters that aren't caught by the 
        universal json applicator (set_GF)
    """
    steering_params.fullLivetime = {0: float(parameters['years'])*365*24*60*60}
    steering_params.sterile_model_label = point
    #steering_params.spline_dom_efficiency = bool(parameters['systematics'][0])
    #steering_params.spline_hole_ice = bool(parameters['systematics'][1])
    #steering_params.load_atmospheric_density_spline = bool(parameters['systematics'][6])
    steering_params.spline_hqdom_efficiency = bool(parameters['systematics'][7])
    #steering_params.load_barr_gradients = '1' in parameters['barr']
    #steering_params.use_ice_gradients = '1' in parameters['multisim']
    steering_params.readCompact = False

steer()

print("------------ > Building GF Fitter")
golemfit = gf.GolemFit(datapaths, steering_params, npp)
print("------------ > Setting Fit Flags")
golemfit.SetFitParametersFlag(fitparams_flag)
golemfit.SetFitParametersPriors(priors)


realization_dist = np.load(run_options['realization'])['realization']
golemfit.Swallow(realization_dist)

# Minimize to the fit
print("------------ > Do Fit")
min_llh = golemfit.MinLLH()
fit_keys = listattr(min_llh.params) # grabs the attributes that we care about (all of them...)

for key in fit_keys:
    try:
        print("{}: {}".format(key, getattr(min_llh.params, key)))
    except AttributeError:
        print("Encountered invalid attribute {} of minimization".format(key))

fit = golemfit.GetExpectation(min_llh.params)[0][0][0] #?????
fit_sum = sum(fit)

output_dict = {}

output_dict['fit_params'] = {}
for key in fit_keys:
    try:
        output_dict['fit_params'][key] = getattr(min_llh.params,key)
    except AttributeError:
        print("Argh. Invalid attribute {}".format(key))

# Write the fit parametrs to a json file
target_file = os.path.join(run_options['outdir'],"GF_fit_" +run_options['point']+ ".json")
with open(target_file,'w') as f:
    json.dump(output_dict, f, ensure_ascii=True, indent=2)
print("-------- > Wrote fit data to {}".format(target_file))
