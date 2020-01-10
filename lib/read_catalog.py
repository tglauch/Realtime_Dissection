import numpy as np

jy_to_erg = 1e-23

def read_from_observation(path):
    name = path.split('/')[-1].split('.')[0]
    if name == 'ovro':
        return read_from_ovro(path) # data format as downloaded from the ovro website
    if name == 'asasn':
        return read_from_asasn(path) #
    if name == 'swift':
        return read_from_swift(path)

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
    flux = idata[:,3]
    print(flux)
    flux_up = (idata[:,3] + idata[:,4])
    flux_low = (idata[:,3] - idata[:,4])
    times = idata[:,0]
    print(times)
    return np.array([frequency, flux, flux_low, flux_up, times]).T
