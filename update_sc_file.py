import time
import wget
from astropy.io import fits
import requests
import os


url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/lat_spacecraft_merged.fits'
current_file = '/scratch8/tglauch/spacecraft_files/current.fits'
new_file = '/scratch8/tglauch/spacecraft_files/new.fits'

def get_sc_file():
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(new_file, 'w+b') as f:
                for chunk in r.iter_content(chunk_size=8192): 
                    if chunk: # filter out keep-alive new chunks
                        f.write(chunk)    
    except Exception as inst:
        print('Faild to Download new Spacecraft data')
        print(inst)
        return False
    if os.path.exists(current_file):
        os.remove(current_file)
    os.rename(new_file, current_file)
    check_file = fits.open(current_file)
    print('\n Photon Data from: {} to {}'.format(check_file[1].header['TSTART'],
                                                 check_file[1].header['TSTOP']))
    check_file.close()
    return True

def check_for_new_sc_file():
    html = requests.get('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/')
    html = html.text.split('lat_spacecraft_merged.fits</a>')[1].split('\n')[0].strip()[:17]
    with open('/scratch8/tglauch/spacecraft_files/last_update.txt', 'r') as f:
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
    
    
