import time
import wget
from astropy.io import fits
import requests
import os 

def get_sc_file():
    try:
        filename = wget.download('ftp://legacy.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/lat_spacecraft_merged.fits',
                                 out='./sc_files/new.fits')
    except Exception as inst:
        print('Faild to Download new Spacecraft data')
        print(inst)
        return False
    new_file = fits.open('./sc_files/new.fits')
    print('\n Photon Data from: {} to {}'.format(new_file[1].header['TSTART'], new_file[1].header['TSTOP']))
    new_file.close()
    os.remove('./sc_files/current.fits')
    os.rename('./sc_files/new.fits', './sc_files/current.fits')
    return True

def check_for_new_sc_file():
    html = requests.get('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/')
    html = html.text.split('lat_spacecraft_merged.fits</a>')[1].split('\n')[0].strip()[:17]
    with open('./sc_files/last_update.txt', 'r') as f:
        last_update = f.read()
    if html != last_update:
        print('UPDATE spacecraft file [modified at {}]'.format(html))
        success = get_sc_file()
        if success == True:
            with open('./sc_files/last_update.txt', 'w') as f:
                f.write(html)
        print('Finished')
    return

while True:
    check_for_new_sc_file()
    print('Wait...')
    time.sleep(30 * 60)
    
    
