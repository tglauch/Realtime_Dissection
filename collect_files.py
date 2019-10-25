import numpy as np
import os
from shutil import copy2
import sys
sys.path.append('./lib')
from lib.read_catalog import read_from_observation 

in_base = '/scratch9/tglauch/realtime_service/output/TXS/'
src = '4FGL_J0509.4+0542'
out_base = '/scratch9/tglauch/erin/'

in_path = os.path.join(in_base, src)
out_path = os.path.join(out_base, src)
if not os.path.exists(out_path):
    os.makedirs(out_path)

c_file = os.path.join(in_path, 'sed.txt')
if os.path.exists(c_file):
    copy2(c_file, out_path)

lc_out = os.path.join(out_path, 'fermi')
if not os.path.exists(lc_out):
    os.makedirs(lc_out)
lc_base = os.path.join(in_path, 'lightcurve')
for f in ['lightcurve.png', 'lightcurve.pdf', 'movie.gif']:
    c_file = os.path.join(lc_base, f)
    copy2(c_file, lc_out)

copy_files = ['bowtie.npy', 'llh.fits', 'sed.fits']
dirs = [i for i in os.listdir(lc_base) if os.path.isdir(os.path.join(lc_base, i))]
for d in dirs:
    if not os.path.exists(os.path.join(lc_out, d)):
        os.makedirs(os.path.join(lc_out, d))
    for f in copy_files:
        copy2(os.path.join(lc_base, d, f),
              os.path.join(lc_out, d))
all_mission_dir = os.path.join(lc_out, 'full_mission')
if not os.path.exists(all_mission_dir):
    os.makedirs(all_mission_dir)
for f in copy_files:
    copy2(os.path.join(in_path, 'all_year/sed', f),
          all_mission_dir)

more_data = os.listdir(os.path.join(in_path, 'add_data'))
for f in more_data:
    add_data = [list(i) for i in read_from_observation(os.path.join(in_path, 'add_data', f))]
    if len(add_data) > 0 :
        for i in  add_data:
            i.extend([i[-1], f.split('.')[0]])
        print add_data
        with open(os.path.join(out_path, 'sed.txt'), 'a') as sed_file:
            np.savetxt(sed_file, add_data, fmt="%s")
