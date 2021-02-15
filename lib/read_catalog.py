import numpy as np

jy_to_erg = 1e-23
keV_to_hz = 2.41799050402417e+17

def read_from_observation(path):
    name = path.split('/')[-1].split('.')[0]
    if  'ovro' in name:
        return read_from_ovro(path) # data format as downloaded from the ovro website
    if 'asasn' in name:
        return read_from_asasn(path) #
    if 'swift' in name:
        return read_from_swift(path)
    if 'xrtproc' in name:
        return read_from_swift_xrtproc(path)
    else:
        return False 

def read_from_ovro(path):
    idata = np.genfromtxt(path, skip_header=1, delimiter=',')
    times = idata[:,0]
    flux = (idata[:,1]) * jy_to_erg * 15 * 1e9
    flux_up = (idata[:,1] + idata[:,2]) * jy_to_erg * 15 * 1e9 
    flux_low = (idata[:,1] - idata[:,2]) * jy_to_erg * 15 * 1e9
    frequency = np.ones(len(flux)) * 15 * 1e9
    return np.array([frequency, flux, flux_low, flux_up, times]).T

def read_from_asasn(path):
    idata = np.genfromtxt(path, delimiter='|')
    times = idata[:,4]
    flux = (idata[:,2])
    flux_up = (idata[:,2] + idata[:,3]) 
    flux_low = (idata[:,2] - idata[:,3])
    frequency = idata[:,0]
    return np.array([frequency, flux, flux_low, flux_up, times]).T

def read_from_swift(path):
    idata = np.genfromtxt(path)
    frequency = idata[:,1]
    flux = idata[:,2]
    flux_up =  idata[:,3]
    flux_low = idata[:,4]
    times = idata[:,0]
    return np.array([frequency, flux, flux_low, flux_up, times]).T


def read_from_swift_xrtproc(path):
    idata = np.genfromtxt(path, delimiter="," , skip_header=1,
                         names=["ra","dec","MJD","F5kev","F5kev_err","F05kev","F05kev_err","F15kev","F15kev_err",
                                "F3kev","F3kev_err","F45kev","F45kev_err","F1kev","F1kev_err","IsDet"],)
    idata = np.atleast_1d(idata)
    frequency = np.array([0.5, 1.5, 4.5]) * keV_to_hz
    flux = np.concatenate([idata["F05kev"], idata["F15kev"], idata["F45kev"]])
    flux_up = np.concatenate([idata["F05kev_err"], idata["F15kev_err"], idata["F45kev_err"]]) + flux
    flux_low = flux - np.concatenate([idata["F05kev_err"], idata["F15kev_err"], idata["F45kev_err"]])
    times = np.concatenate([[i]*3 for i in idata['MJD']]) 
    frequency = np.concatenate([[i]*len(idata) for i in [0.5, 1.5, 4.5]]) * keV_to_hz
    print(np.array([frequency, flux, flux_low, flux_up, times]).T)
    return np.array([frequency, flux, flux_low, flux_up, times]).T
