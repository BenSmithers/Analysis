import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

import h5py as h5

# load data
data_raw = h5.File("simple_data.hdf5",'r')
NC_data = np.array(data_raw['NC'])
CC_data = np.array(data_raw['CC'])
data_raw.close()

NC_data = np.transpose(NC_data)
CC_data = np.transpose(CC_data)

energies = np.logspace(2, 8, 21)
plt.hist(NC_data[0],bins=energies, weights=NC_data[1], alpha=0.5, label="NC")
plt.hist(CC_data[0],bins=energies, weights=CC_data[1], alpha=0.5,label="CC")
plt.xscale('log')
plt.legend()
plt.yscale('log')
plt.xlabel("Energy [GeV]",size=14)
plt.savefig("test.png",dpi=400)
