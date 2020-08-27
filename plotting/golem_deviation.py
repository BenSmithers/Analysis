import numpy as np
import json
import os
from math import inf

import matplotlib
import matplotlib.pyplot as plt

#datadir = "/data/user/bsmithers/runs"
datadir = "/home/benito/software/data"

#bpl_raw = "GF_fit_BPL_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.json"
bpl_raw =  "BPL0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.json"
#pl_raw = "GF_fit_PL_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.json"
pl_raw =  "PL0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.json"

# load the files in
if not os.path.exists(os.path.join(datadir, bpl_raw)):
    raise IOError("Broken Power Law file not found: {}".format(bpl_raw))

if not os.path.exists(os.path.join(datadir,pl_raw)):
    raise IOError("Power Law file not found {}".format(pl_raw))

fl = open(os.path.join(datadir, bpl_raw), 'r')
bpl_data = json.load(fl)
fl.close()

fl = open(os.path.join(datadir, pl_raw),'r')
pl_data = json.load(fl)
fl.close()

# make a little dictionary to make this cleaner
data = {"Broken Power Law": bpl_data, "Single Power Law": pl_data}

# using this to keep track of the heights AND the values 
data_heights = {}

current_x = 0
data_labels = {}

for datakey in data.keys():
    data_heights[datakey]={}
    label_done = False

    for key in data[datakey]["fit_params"].keys():
        assert(isinstance(key, str))

        if key+"Center" not in data[datakey]["priors"]:
            print("Skipping '{}'".format(key))
            continue
        if (data[datakey]["fit_params"][key] - data[datakey]["priors"][key+"Center"])!=0:        
            data_heights[datakey][key] = data[datakey]["fit_params"][key] - data[datakey]["priors"][key+"Center"]
            
            if True: #data[datakey]["priors"][key+"Width"] != inf:
                data_heights[datakey][key]/= data[datakey]["priors"][key+"Width"]
            else:
                print("Not dividing that inf!")
        else:
            continue

        if key not in data_labels:
            data_labels[key] = current_x
            current_x += 1

        if datakey=="Broken Power Law":
            color = 'r'
        else:
            color = 'b'

        plt.plot(data_heights[datakey][key], data_labels[key]  ,marker="o", ms=10, ls='', color=color,alpha=0.5, label=(None if label_done else datakey))
        label_done = True

plt.yticks( list(data_labels.values()), list(data_labels.keys())) #, rotation=90)
plt.vlines(0,ymin=-0.5, ymax=current_x-0.5, alpha=0.5, zorder=-1)
plt.title("Pull / Width")
plt.ylim([-0.5, current_x-0.5])
plt.tight_layout()
plt.legend()
plt.grid(which='both',axis='y',alpha=0.7)
plt.savefig("deviation.png",dpi=400)
plt.show()
