## import GolemFitPy as gf # central fit utility

import json # used to load in a defaults file 
import os 

import numpy as np

configuration_file = "configuration.json"
# using os to open the file correctly, regardless of how this script is called
f = open(os.path.join(os.path.dirname(__file__), configuration_file), 'r')
config = json.load(f)
f.close()

run_options= config['run_options']
fit_params = config['central_values']
parameters = config['parameters']

def parse_point( this_point ):
    if not isinstance(this_point, str):
        raise TypeError("Expected {}, not {}".format(str, type(this_point)))
    broken = this_point.split("_")
    if not (len(broken)==7):
        raise ValueError("{} appears not to be a valid point. Expected 7 parts.".format(this_point))
    
    def round_string(number):
        return("{:.4f}".format(np.round(float(number),4)))

    dm2 = round_string(broken[1])
    th24= float(broken[3])
    s2t = "{:.4f}".format(np.sin(2*th24)**2)
    th34= round_string(broken[4])
    Uu4sq= np.sin(th24)**2
    Ut4sq= (np.cos(th24)**2)*(np.sin(float(th34))**2)
    return( Uu4sq, Ut4sq )


# transcribe the fit parameters from the json file into this object 
fitparams = gf.FitParameters(gf.sampleTag.Sterile)
for param in fit_params.keys():
    setattr(param, fitparams, fit_params[param]) 

steering_params = gf.SteeringParams(gf.sampleTag.Sterile)
