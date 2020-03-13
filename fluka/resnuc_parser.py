from glob import glob # load files 

from flutils import *

from numpy import mean, std, sqrt, log10
import matplotlib.pyplot as plt
import numpy as np

#print("Data File Parsed: {}".format(list_shape( data_array)))


all_files = glob( "/home/benito/software/fluka/runs/xenon-136/*_fort*" )

volume = 7.8540e5
maxz = 60
maxn_z = 85
minn_z = -4
k = -5

results = {}

def get_elm( Z ):
    elms = {52:"Te", 53:"I", 54:"Xe", 55:"Cs"}
    if not isinstance(Z, int):
        raise TypeError("Should be {}, got {}".format(int, type(Z)))
    
    try:
        return( elms[Z] )
    except KeyError:
        return( "??" )
    

def add( key, value):
    if key not in results:
        results[key] = [ value ]
    else:
        results[key] += [ value ]


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
