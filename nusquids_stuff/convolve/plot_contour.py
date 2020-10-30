import pickle 
import os 

from utils import get_loc, bhist, bilinear_interp

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt 

savefile_raw = ".analysis_level.dat"
savefile = os.path.join(os.path.dirname(__file__), savefile_raw)

def load_file():
    if not os.path.exists(savefile):
        raise Exception("Datafile not found. Generate data before calling this")
    print("Loading Data")

    f = open(savefile,'rb')
    all_data = pickle.load(f)
    f.close()
    
    e_reco = all_data["e_reco"]
    e_true = all_data["e_true"]
    a_reco = all_data["a_reco"]
    a_true = all_data["a_true"]
    probs = all_data["flux"]

    return(e_reco, e_true, a_reco, a_true, probs)

def build_contours(obs_energy, obs_angle):
    """
    Presuming we're given some observed angle and observed energy...

    We first figure out which event bin we want 
    """

    # these are bin edges
    e_reco, e_true, a_reco, a_true, probs = load_file()

    # let's get the centers we want
    e_true_centers = bhist([e_true]).centers
    a_true_centers = bhist([a_true]).centers
    
    e_reco_centers = bhist([e_reco]).centers
    a_reco_centers = bhist([a_reco]).centers
    
    e_left, e_right = get_loc(obs_energy, e_reco_centers)
    a_left, a_right = get_loc(obs_angle, a_reco_centers)

    # so this is a bit wild. The reco-space entries are coarse. So, we grab neighboring 2D reconstruction arrays and do a bilinear interpolation of 2D arrays around our given point in reco-space 
    key = list(probs.keys())[0]
    p0=(obs_energy, obs_angle)
    p1=(e_reco_centers[e_left], a_reco_centers[a_left])
    p2=(e_reco_centers[e_right],a_reco_centers[a_right])
    q11=probs[key][e_left][:][a_left][:]
    q21=probs[key][e_right][:][a_left][:]
    q12=probs[key][e_left][:][a_right][:]
    q22=probs[key][e_right][:][a_right][:]

    focus_prob = bilinear_interp(p0,p1,p2,q11,q12,q21,q22)
    del probs # clear this out of memory. It's really big... 

    # pcolormesh 
    
    plt.pcolormesh( e_true_centers, a_true_centers, focus_prob)
    plt.xlabel("True Energy [GeV]",size=14)
    plt.xscale('log')
    plt.ylabel("Cos Zenith", size=14)
    plt.title(r"Expected Truth for E={}, $\cos\theta$={}".format(obs_energy, obs_angle))
    plt.savefig("contour.png",dpi=400)
    plt.show()

build_contours(10**13, -0.9)
