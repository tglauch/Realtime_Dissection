import numpy as np
import os
from shutil import copy2, make_archive
import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('./lib')
from lib.read_catalog import read_from_observation 
import zipfile
import argparse
import pickle
from analysis import * 

def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file))

p = argparse.ArgumentParser(
       description = "collecting realtime pipeline output",
       formatter_class = argparse.RawTextHelpFormatter)
p.add_argument("--bpath", type=str,
     help='What is the basepath')
p.add_argument("--src", type=str, required=True,
     help='Which source?')
args = p.parse_args()
print('Run with args')
print(args)

in_base = os.path.normpath(args.bpath)
src = args.src
analysis_object_path = os.path.join(in_base,'analysis.pickle')
with open(analysis_object_path, "rb") as f:
    analysis = pickle.load(f)
    mjd = analysis.mjd
add_name = in_base.split('/')[-1]
out_base = os.path.join('/scratch8/tglauch/all_blazar_dissection_new3', add_name)

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

copy_files = ['bowtie.npy', 'llh.fits', 'sed.fits', 'llh.npy', 'sed.npy', 'sed.pdf']
dirs = [i for i in os.listdir(lc_base) if os.path.isdir(os.path.join(lc_base, i))]
for d in dirs:
    if not os.path.exists(os.path.join(lc_out, d)):
        os.makedirs(os.path.join(lc_out, d))
    for f in copy_files:
        try:
            copy2(os.path.join(lc_base, d, f),
                  os.path.join(lc_out, d))
        except Exception as e:
            print('{} does not exist'.format(os.path.join(lc_base, d, f)))
    mjds = d.split('_')
    if (float(mjds[0]) < mjd) & (float(mjds[1]) > mjd):
        copy2(os.path.join(lc_base, d, 'sed.pdf'),
              os.path.join(out_path))
all_mission_dir = os.path.join(lc_out, 'full_mission')
if not os.path.exists(all_mission_dir):
    os.makedirs(all_mission_dir)
for f in copy_files:
    copy2(os.path.join(in_path, 'all_year/sed', f),
          all_mission_dir)

more_data = os.listdir(os.path.join(in_path, 'add_data'))
#for f in more_data:
#    add_data = [list(i) for i in read_from_observation(os.path.join(in_path, 'add_data', f))]
#    if len(add_data) > 0 :
#        for i in  add_data:
#            i.extend([i[-1], f.split('.')[0]])
#        with open(os.path.join(out_path, 'sed.txt'), 'a') as sed_file:
#            np.savetxt(sed_file, add_data, fmt="%s")
zip_path = os.path.normpath(out_path)+'.zip'
make_archive(out_path, 'zip', out_path)

#print('Create Zip File in {}'.format(zip_path))
#zipf = zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED)
#zipdir(out_path, zipf)
#zipf.close()
