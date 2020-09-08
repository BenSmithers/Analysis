import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import os
from math import sqrt

emin = 1
emax = 150
n_bins = 50

def get_invmass_from_data( data ):

    # invariant mass is 
    # sqrt( m1^2 + m2^2 + 2(E1*E2 + p1\dot p2) )
    m1sq = np.array([ evt[5]**2 - evt[2]**2 - evt[3]**2 - evt[4]**2  for evt in data])
    m2sq = np.array([ evt[11]**2 -evt[10]**2 -evt[9]**2 - evt[8]**2  for evt in data])
    p1dp2 = 2*np.array([ evt[11]*evt[5] - evt[2]*evt[8] - evt[3]*evt[9] - evt[4]*evt[10] for evt in data])

    invmass = np.sqrt(m1sq + m2sq + p1dp2)

    return(invmass)

def plot_files(files):

    for each in files:
        if not os.path.exists(each):
            print("File '{}' does not exist, skipping".format(each))
            continue

        # load the event
        data = np.loadtxt(each)
        invmass = get_invmass_from_data(data)

        bins = np.logspace(np.log10(emin), np.log10(emax), n_bins+1)
#        bins = np.linspace(emin, emax, n_bins+1)

        plt.hist(invmass, bins,alpha=0.5, label=each.split(".")[0])

    plt.legend()
    plt.title("Ben's Plot - Do Not Steal")
    plt.xlabel("Invariant Mass [GeV]",size=14)
    plt.ylabel("Counts",size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    plt.close()

# load in flavor pairs
files_flavor = ["mu_el.txt", "two_el.txt", "two_mu.txt"]

files_sign = ["same_sign.txt", "opp_sign.txt"]

plot_files(files_flavor)
plot_files(files_sign)

