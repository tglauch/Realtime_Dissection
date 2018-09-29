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
vou_out = os.path.join(bpath, 'vou_blazar' )
if not os.path.exists(vou_out):
    os.makedirs(vou_out)
os.chdir(vou_out)
subprocess.call(cmd)
print('Get Sources....')
src_dict = get_sources(args['ra'], args['dec'])
os.chdir(this_path)
get_data(args['ra'], args['dec'], emin=2000,
         dt=364, out_dir=fermi_data)  # dt hardcoded!!!!
for i, src in enumerate(src_dict['name']):
    bpath_src = os.path.join(bpath, src.replace(' ', '_'))
    print bpath_src
    if os.path.exists(bpath_src):
        shutil.rmtree(bpath_src)
    os.makedirs(bpath_src)
    args = '--target_src {} --free_radius 2 --data_path {} --use_3FGL --outfolder {} '.format(src.replace(' ', '_'), fermi_data, bpath_src)
    if '3FGL' not in src:
        xml_path = os.path.join(bpath_src, 'add_source.xml')
        generate_src_xml(src, src_dict['ra'][i], src_dict['dec'][i], xml_path)
        args += '--xml_path {}'.format(xml_path)

    with open('./slurm_draft.xml', 'r') as f:
        submit = (f.read()).format(bpath=bpath_src, args=args)
    submitfile = os.path.join(bpath_src, 'submit.sub')
    with open(submitfile, "w+") as file:
        file.write(submit)
    os.system("sbatch {}".format(submitfile))
