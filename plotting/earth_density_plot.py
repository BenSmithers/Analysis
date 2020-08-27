import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import numpy as np

import os
param_dir = "/home/benito/software/LeptonInjector/source/resources/earthparams/densities"
prem_file = "PREM_mmc.dat"

datafile = os.path.join(param_dir, prem_file)

def parse_line(line):
    """
    Parses the line in the data file, returns a function that returns the density at a given radius 
    """
    pass

def trim_whitespace(line):
    """
    Gets rid of any leading whitespace, drops out comments and stuff 
    """
    if not isinstance(line, str):
        raise TypeError("Expected {}, got {}".format(str, type(line)))

    temp = line[:]
    while temp[0]==" ":
        temp = temp[1:]

    return(temp)

def clear_comments(line):
    """
    scans for comment characters...
    """

    comment = "#"
    iterate = 0
    while iterate<len(line):
        if line[iterate]==comment:
            return(line[:iterate])

        iterate+=1
    return(line)

data = open(datafile, 'r')

parsed_lines = []
for line in data:
    trimmed = clear_comments(line)
    if trimmed=="":
        continue
    trimmed = trim_whitespace(trimmed)

    if trimmed[-1]=="\n":
        trimmed = trimmed[:-1]
    if trimmed=="":
        continue

    parsed_lines.append(trimmed.split(None))


def build_density_function(*args): 
    def density(radius):
        if len(args)==1:
            return((radius**0)*float(args[0]))
        else:
            value = 0
            for term in range(len(args)):
                value += (radius**term)*float(args[term])
            return(value)

    return(density)


color = "hot"
this_cmap = plt.get_cmap(color)
def get_color(n, colormax=6.0):
    #0,1,2,3,
    colorVal = this_cmap(1.-n/colormax)
    return(colorVal) # (r, g, b)


plt.figure(1)
last_rad = 0.
lineno = 0

active_name = ""
active_index = 0
for line in parsed_lines:
    new_legend = False
    outer_rad = float(line[0])
    name_1 = line[1]
    name_2 = line[2]
    some_param = line[3]

    if name_2!=active_name:
        active_name = name_2
        active_index+=1
        new_legend = True

    rest = tuple(line[4:])

    density_func = build_density_function(*rest)

    sample_radii = np.logspace(np.log10(6478001 - outer_rad), np.log10(6478001-last_rad), 100)
    #sample_radii = np.logspace(np.log10(last_rad), np.log10(outer_rad), 100)
    sample_densities = density_func(6478001 - sample_radii)

    last_rad = outer_rad

    plt.plot(sample_radii, sample_densities, color=get_color(active_index), label=(active_name if new_legend else None), lw=4)
    lineno+=1

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Depth [cm]",size=14)
plt.ylabel("Density [g/cm3]",size=14)
plt.ylim([0.0001, 20])
plt.xlim([100000, 6500000])
plt.legend(loc=4)
plt.show()
