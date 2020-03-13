from glob import glob # load files 

from flutils import *

from numpy import mean, std, sqrt, log10
import matplotlib.pyplot as plt
import numpy as np

"""
Ben Smithers
benjamin.smithers@mavs.uta.edu

This script loads a group of outputs from the RESNUCLEi card in a fluka run and plots the residual nuclei on a bar graph. 

Users need to specify values from the runs' headers. Assumes that these headers are all the same.

TODO: allow script to automatiaclly read the header! 
"""

#print("Data File Parsed: {}".format(list_shape( data_array)))


all_files = glob( "/home/benito/software/fluka/runs/xenon-136/*_fort*" )

volume = 7.8540e5
maxz = 60
maxn_z = 85
minn_z = -4
k = -5

# We score all results in this dictionary. 
# Each key is an isotope ("Xe136"), and the value is the number of nuclei/primary 
# the value is a list of such quantities: one value for each registered run
results = {}

def get_elm( Z ):
    """
    Returns the 1-2 letter combo corresponding to a element. Obviously could just be a dictionary, but I wanted to handle the exceptions. 

    Todo: look into database/library already existing with this dictionary's information 
    """
    elms = {52:"Te", 53:"I", 54:"Xe", 55:"Cs", 1:"H", 2:"He"}
    if not isinstance(Z, int):
        raise TypeError("Should be {}, got {}".format(int, type(Z)))
    
    try:
        return( elms[Z] )
    except KeyError:
        print("Add Z={}!".format(Z))
        return( "??" )
    

def add( key, value):
    """
    Registers the given flux with the given key (isotope) 
    """
    global results
    if key not in results:
        results[key] = [ value ]
    else:
        results[key] += [ value ]

# load in each data file, parse it, and register the entries with results dict 
for each in all_files:
    data = load_datafile( each )
    reform = reshape( data, [ maxn_z - minn_z +1, maxz ] )


    iter_z = 1
    while iter_z < maxz:
        iter_what = minn_z - k
        while iter_what < (maxn_z-k):
            A = iter_what + k +2*iter_z
            numb = reform[iter_what-1][iter_z-1]
            if numb!=0:
                add( "{}_{}".format( iter_z, A), numb*volume )
#                print("{} {} {}".format(iter_z, A, numb*volume ))
            iter_what += 1

        iter_z+=1 

values = []
errors = []
labels = []
for key in results:
    desc = key.split("_")

    element = get_elm( int(desc[0]) )

    order = int(log10( mean(results[key]) ))-1
    values.append( mean(results[key]) )
    errors.append( std(results[key])/sqrt(len(results[key])) ) 
    labels.append( element+str(desc[1]) )

    print( element + str(desc[1]) + r": ({:.2f} +/- {:.2f}) x 10^{}".format(mean(results[key])/(10**order), std(results[key])/(sqrt(len(results[key]))*(10**order) ) , order ) )

positions = range( len(values))
fig, ax = plt.subplots()
ax.bar(x=positions, height=values, width=0.8, yerr=errors, align='edge', ecolor='black',capsize=10)
ax.set_xticks( np.array(positions)+0.4 )
ax.set_xticklabels( labels, size=14 )
ax.set_ylabel("Isotopes per Neutron", size=16)
ax.set_yscale('log')
ax.text( 0, max(values)*0.9 , "FLUKA Preliminary", fontsize=16, color=(1.0,0,0))
plt.show()
