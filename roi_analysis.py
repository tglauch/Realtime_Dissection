#! /usr/bin/env python
# coding: utf-8

import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('/scratch9/tglauch/realtime_service/python/lib/python2.7/site-packages')
from get_fermi_data import get_data
import subprocess
import argparse
import os
from roi_functions import get_sources
import shutil
import datetime
from slack_lib import print_to_slack
import numpy as np

MJDREF = 51910.0 + 7.428703703703703E-4
def MET_to_MJD(met_time):
    return float(met_time / 86400. + MJDREF)

def parseArguments():
    """Parse the command line arguments
    Returns:
    args : Dictionary containing the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ra",
        help="right ascension",
        type=float, required=True)
    parser.add_argument(
        "--dec",
        help="declination",
        type=float, required=True)
    args = parser.parse_args()
    return args.__dict__


def submit_fit(args, bpath, src_dict, which, trange=''):
    if trange is '':
        opath = os.path.join(bpath, 'all_year')
        args += ' --outfolder {} '.format(opath)
    else:
        opath = os.path.join(
            bpath, 'lightcurve/{:.1f}_{:.1f}'.format(trange[0], trange[1]))
        args += ' --time_range {} {} '.format(trange[0], trange[1])
        args += ' --outfolder {} '.format(opath)
    os.makedirs(opath)
    if '3FGL' not in src:
        xml_path = os.path.join(opath, 'add_source.xml')
        generate_src_xml(src_dict['name'][which], src_dict['ra'][which],
                         src_dict['dec'][which], xml_path)
        args += '--xml_path {}'.format(xml_path)
    with open('./slurm_draft.xml', 'r') as f:
        submit = (f.read()).format(bpath=opath, args=args)
    submitfile = os.path.join(opath, 'submit.sub')
    with open(submitfile, "w+") as file:
        file.write(submit)
    os.system("sbatch {}".format(submitfile))


def generate_src_xml(name, ra, dec, xml_path):
    with open('./src_xml_draft.xml', 'r') as f:
        xml_str = (f.read()).format(name=name, ra=ra, dec=dec)
    with open(xml_path, 'w+') as f:
        f.write(xml_str)
    return


dtime = datetime.datetime.now()
ev_str = 'IC{}{:02d}{}'.format(str(dtime.year)[-2:], dtime.month, dtime.day)
bpath = '/scratch9/tglauch/realtime_service/output/{}'.format(ev_str)
this_path = os.path.dirname(os.path.abspath(__file__))
fermi_data = os.path.join(bpath, 'fermi_data')


args = parseArguments()
cmd = [os.path.realpath('/scratch9/tglauch/VOU_Blazars/bin/vou-blazars'),
       str(args['ra']), str(args['dec']), str(120)]
vou_out = os.path.join(bpath, 'vou_blazar')
if not os.path.exists(vou_out):
    os.makedirs(vou_out)
os.chdir(vou_out)
subprocess.call(cmd)
print('Get Sources....')
src_dict, out_str = get_sources(args['ra'], args['dec'])
cand_ps = os.path.join(vou_out, 'candidates.ps')
cand_png= os.path.join(vou_out, 'candidates.png')
os.system('convert ' + cand_ps + ' -density 600 ' + cand_png )
print_to_slack(out_str, cand_png)
'''
os.chdir(this_path)
MET = get_data(args['ra'], args['dec'], emin=2000,
               dt=56, out_dir=fermi_data)  # dt hardcoded!!!!
MJD = [MET_to_MJD(float(i)) for i in MET]
print MJD
for i, src in enumerate(src_dict['name']):
    bpath_src = os.path.join(bpath, src.replace(' ', '_'))
    print bpath_src
    if os.path.exists(bpath_src):
        shutil.rmtree(bpath_src)
    args = '--target_src {} --free_radius 2 --data_path {} --use_3FGL '
    print fermi_data
    args = args.format(src.replace(' ', '_'), fermi_data)
    print args
    time_windows = [[k, k+28.] for k in np.arange(MJD[0], MJD[1], 28)]
    time_windows.append('')
    for t_window in time_windows:
        print t_window
        submit_fit(args, bpath_src, src_dict, i, trange=t_window)

'''
