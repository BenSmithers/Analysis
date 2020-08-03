import LeptonWeighter as LW

from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube import LeptonInjector 

import os
from glob import glob
import h5py

root_dir = "/data/user/eganster/GlobalFit/simulation/datasets/SnowStorm_NuTau/005"

infiles = glob(os.path.join(root_dir,"gen_*.i3.zst"))
ignore = "poly"
infiles = list(filter(lambda x: ignore not in x, infiles))


print("Loading: {}".format(infiles))
lic_files = glob(os.path.join(root_dir,"gen_*.lic"))
print("Found LIC Files: {}".format(lic_files))


net_gen = []
for lic in lic_files:
    net_gen += LW.MakeGeneratorsFromLICFile( lic )

flux_params={ 'constant': 10**-18, 'index':-2, 'scale':10**5 }
liveTime   =3.1536e7
flux       = LW.PowerLawFlux(flux_params['constant'] , flux_params['index'] , flux_params['scale'] )
xs         = LW.CrossSectionFromSpline(
            "/data/user/bsmithers/cross_sections/dsdxdy_nu_CC_iso.fits",
            "/data/user/bsmithers/cross_sections/dsdxdy_nubar_CC_iso.fits",
            "/data/user/bsmithers/cross_sections/dsdxdy_nu_NC_iso.fits",
            "/data/user/bsmithers/cross_sections/dsdxdy_nubar_NC_iso.fits")
weight_event = LW.Weighter(flux,xs,net_gen)

def get_weight(frame):
    LWevent = LW.Event()
    EventProperties                 = frame['EventProperties']

    # get MCPrimary
    if "MCPrimary1" in frame:
        MCPrimary                   = frame["MCPrimary1"]
    else:
        MCPrimary                   = frame["I3MCTree"].primaries[0]
        frame["MCPrimary1"] = MCPrimary

    LeptonInjectorProperties        = frame['LeptonInjectorProperties']
    LWevent.primary_type            = LW.ParticleType(EventProperties.initialType)
    LWevent.final_state_particle_0  = LW.ParticleType(EventProperties.finalType1)
    LWevent.final_state_particle_1  = LW.ParticleType(EventProperties.finalType2)
    LWevent.zenith                  = EventProperties.zenith
    LWevent.azimuth                 = EventProperties.azimuth
    LWevent.energy                  = EventProperties.totalEnergy
    LWevent.interaction_x           = EventProperties.finalStateX
    LWevent.interaction_y           = EventProperties.finalStateY

    LWevent.x                       = MCPrimary.pos.x
    LWevent.y                       = MCPrimary.pos.y
    LWevent.z                       = MCPrimary.pos.z

    if type(EventProperties)==LeptonInjector.VolumeEventProperties:
        LWevent.total_column_depth = EventProperties.totalColumnDepth
        LWevent.radius = EventProperties.radius
    else:
        LWevent.total_column_depth  = EventProperties.totalColumnDepth
        LWevent.radius              = EventProperties.impactParameter


    weight_OneWeight                = weight_event.get_oneweight(LWevent)
    weight_Powerlaw                 = weight_event(LWevent)
    
    return(weight_OneWeight)

outdata_NC = []
outdata_CC = []

for infile in infiles:
    write_CC = "CC" in infile    

    data = dataio.I3File(infile, 'r')
    while data.more():
        frame = data.pop_frame()
        if str(frame.Stop)!='DAQ':
            continue
        if write_CC:
            outdata_CC.append((frame['EventProperties'].totalEnergy,get_weight(frame)))
        else:
            outdata_NC.append((frame['EventProperties'].totalEnergy,get_weight(frame)))

write_dir = os.path.dirname(__file__)
file_name = os.path.join(write_dir, "simple_data.hdf5")
outdict = {'NC':outdata_NC, "CC":outdata_CC}

phil = h5py.File(file_name, 'w')
for key in outdict:
    dset = phil.create_dataset(key, data=outdict[key])
phil.close()
