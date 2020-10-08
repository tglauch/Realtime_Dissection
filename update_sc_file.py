import time
import wget
from astropy.io import fits
import requests
import os


url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/lat_spacecraft_merged.fits'
current_file = '/scratch8/tglauch/spacecraft_files/current.fits'
new_file = '/scratch8/tglauch/spacecraft_files/new.fits'
last_update_file = '/scratch8/tglauch/spacecraft_files/last_update.txt'


def get_sc_file():
    if os.path.exists(new_file):
        os.remove(new_file)
    try:
        response = requests.head(url)
        print(response.headers)
        c_length = response.headers['Content-Length']
        print('Download Size: {} Bytes'.format(c_length))
        with requests.get(url, stream=True, timeout=30) as r:
            r.raise_for_status()
            with open(new_file, 'w+b') as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024 * 750):
                    f.write(chunk)
    except Exception as inst:
        print('Failed to download new spacecraft data')
        print(inst)
        return False
    dl_size = os.path.getsize(new_file) 
    print('Downloaded {} Bytes'.format(dl_size))
    if (r.status_code == 200) & (int(dl_size) == int(c_length)):
        print('Download Successfull')
    else:
        print('Download Failed. Status: {}; {} ?= {}'.format(r.status_code, dl_size, c_length))
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
    try:
        html = requests.get('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/')
    except Exception as inst:
        print('Failed to get spacecraft info data')
        print(inst)
        return False
    html = html.text.split('lat_spacecraft_merged.fits</a>')[1].split('\n')[0].strip()[:17]
    with open(last_update_file, 'r') as f:
        last_update = f.read()
    if html != last_update:
        print('UPDATE spacecraft file [modified at {}]'.format(html))
        success = get_sc_file()
        if success == True:
            with open(last_update_file, 'w') as f:
                f.write(html)
        else:
            print('File not updated...try again')
            check_for_new_sc_file()
        print('Finished')
    return

while True:
    check_for_new_sc_file()
    print('Wait...')
    time.sleep(60 * 60)
    
    
