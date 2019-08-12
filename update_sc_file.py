import time
import wget
from astropy.io import fits
import requests

def get_sc_file():
    time.sleep(1)
    filename = wget.download('ftp://legacy.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/lat_spacecraft_merged.fits',
                             out='./sc_files/new.fits')
    new_file = fits.open('./sc_files/new.fits')
    print('Photon Data from: {} to {}'.format(x[1].header['TSTART'], x[1].header['TSTOP']))
    os.remove('/sc_files/current.fits')
    os.rename('./sc_files/new.fits', './sc_files/current.fits')
    return

def check_for_new_sc_file():
    html = requests.get('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/')
    html = html.text.split('lat_spacecraft_merged.fits</a>')[1].split('\n')[0].strip()[:17]
    with open('./sc_files/last_update.txt', 'r') as f:
        last_update = f.read()
    if html != last_update:
        print('UPDATE spacecraft file [modified at {}]'.format(html))
        get_sc_file()
        with open('./sc_files/last_update.txt') as f:
            f.write(html)
    return

while True:
    check_for_new_sc_file()
    time.sleep(30 * 60)
    
    
