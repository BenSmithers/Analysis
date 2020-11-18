import nuSQUIDSpy as nsq

# define a couple utility functions
def get_flavor( key ):
    '''
    take a flux dictionary key and return the nusquids flavor type
    
    The dictionary key will be like "electorn_stuff_stuff"
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))

    part = key.split('_')[0].lower()
    if part in ['e', 'eleectron']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.electron )
    elif part in ['mu', 'muon']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.muon )
    elif part in ['tau']:
        return( nsq.NeutrinoCrossSections_NeutrinoFlavor.tau )
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))

def get_neut( key ):
    '''
    Takes a flux dictionary key and returns the nusquids neutrino type 
    (anti neutrino or vanilla neutrino)
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))   
    part = key.split('_')[1].lower()

    if part in ['nu', 'matter','neutrino']:
        return(nsq.NeutrinoCrossSections_NeutrinoType.neutrino)
    elif part in ['nubar', 'antimatter', 'antineutrino']:
        return(nsq.NeutrinoCrossSections_NeutrinoType.antineutrino)
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))

def get_curr(key):
    '''
    Takes a flux dictionary key and returns the nusquids neutrino type 
    '''
    if not isinstance(key, str):
        raise TypeError("Expected {}, got {}".format(str, type(key)))   
    part = key.split('_')[2].lower()

    if part in ['neutral', 'nc']:
        return(nsq.NeutrinoCrossSections_Current.NC)
    elif part in ['charged', 'cc']:
        return(nsq.NeutrinoCrossSections_Current.CC)
    else:
        raise ValueError("Not sure how to work with {}, extracted from {}".format(part, key))


