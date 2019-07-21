# coding: utf-8

import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('/scratch9/tglauch/realtime_service/python/lib/python2.7/site-packages')
sys.path.append('/home/ga53lag/Software/python_scripts/')
import subprocess
import argparse
import os
from roi_functions import submit_fit, partition_t, get_68_psf 
import shutil
import datetime
import numpy as np
import plot
import warnings
import re
import time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
from collections import OrderedDict
from analysis import Analysis
import pickle
from slack_lib import print_to_slack


def parseArguments():
    """Parse the command line arguments
    Returns:
    args : Dictionary containing the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ra",
        help="right ascension",
        type=float,)
    parser.add_argument(
        "--dec",
        help="declination",
        type=float,)
    parser.add_argument(
        "--mjd",
        help="The MJD time of the event",
        type=float)
    parser.add_argument(
        "--mjd_range",
        help="The mjd range for the analysis",
        type=float, nargs='+')
    parser.add_argument(
        "--dt_lc",
        help="time lenght of bins in the light curve",
        type=float, default=200)
    parser.add_argument(
        "--dt",
        help="length of time window to analyze",
        type=float)
    parser.add_argument(
        "--radius",
        help="radius of the region to analyze",
        type=str, default="180")
    parser.add_argument(
        "--emin",
        help="lower energy bound for SED",
        type=float, default=1000)
    parser.add_argument(
        "--event",
        help="event name",
        default='None')
    parser.add_argument(
        "--only_vou",
        help="event name",
        action="store_true", default = False)
    parser.add_argument(
        "--gcn",
        help="pass a link to a gcn notice",
        type=str)
    parser.add_argument(
        "--mode",
        help="end if given mjd is at the end of the time window, mid if in the middle'",
        type=str, default='end')
    parser.add_argument(
        "--recovery",
        help="Is this a run that recovers the previous processing?",
        action="store_true", default = False)
    parser.add_argument(
        "--make_pdf",
        help="Only create the output pdf form a previous run ",
        action="store_true", default = False)
    parser.add_argument(
        "--overwrite",
        help="Only create the output pdf form a previous run ",
        action="store_true", default = False)
    parser.add_argument(
        "--basepath",
        help = 'basepath for the output',
        type=str, default = '/scratch9/tglauch/realtime_service/output/')
    parser.add_argument(
        "--max_dist",
        help = 'maximum distance for the fermi analysis',
        type=float, default=2.) 
    args = parser.parse_args()
    return args.__dict__

def src_path(bpath, src):
    return os.path.join(bpath, src.replace(' ', '_'))


args = parseArguments()
make_pdf = args['make_pdf']
print('Run with args \n {}'.format(args))

# Setup Analysis Class
bpath = os.path.join(args['basepath'], args['event'])
analysis_object_path = os.path.join(bpath,'analysis.pickle')
if os.path.exists(analysis_object_path):
    print('Folder {} already exist....exit. Only create PDF'.format(bpath))
    with open(analysis_object_path, "rb") as f:
        analysis = pickle.load(f)
    analysis.update_gcn()
    make_pdf=True
else:
    ## read GCN
    if args['gcn'] is not None:
        analysis = Analysis()
        analysis.from_gcn(args['gcn'])
    else:
        if args['mjd'] is not None:
            t = Time(args['mjd'], format='mjd')
        else:
            t = Time.now()
            args['mjd'] = t.mjd
        analysis = Analysis(args['mjd'], args['ra'], args['dec'])

    if args['event'] is not 'None':
        ev_str = args['event']
    else:
        t = Time(analysis.mjd, format='mjd')
        ev_str = 'IC{}{:02d}{:02d}'.format(str(t.datetime.year)[-2:],
                                           t.datetime.month, t.datetime.day)
    bpath = os.path.join(args['basepath'], ev_str)
    analysis_object_path = os.path.join(bpath,'analysis.pickle')
    analysis.bpath = bpath
    analysis.emin = args['emin']
    analysis.event_name = ev_str
    if 'err90' in args.keys():
        analysis.err90 = args['err90']
this_path = os.path.dirname(os.path.abspath(__file__))
analysis.this_path = this_path

if not args['recovery'] and not make_pdf:
    # Run VOU Tool
    analysis.ROI_analysis(args['radius'])
    if args['only_vou']:
        exit()

    # Download gamma-ray data
    analysis.get_fermi_data(days=args['dt'], mjd_range=args['mjd_range'])
    with open(analysis_object_path, "wb") as f:
        pickle.dump(analysis, f)



# Start the gamma-ray analysis

if not make_pdf:
    args['overwrite'] = True
    if args['mode'] == 'end':
        tsmjd1 = analysis.mjd-200
        tsmjd2 = analysis.mjd

    else:
        tsmjd1 = analysis.mjd-100
        tsmjd2 = analysis.mjd+100
    analysis.tsmjd1 = tsmjd1
    analysis.tsmjd2 = tsmjd2
    ts_emin = np.max([1000, analysis.emin])
    analysis.ts_emin = ts_emin
    # TS maps
    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
    sargs = sargs.format(get_68_psf(ts_emin), analysis.fermi_data, analysis.ts_emin, analysis.ra, analysis.dec)
    ts_map_path = os.path.join(analysis.bpath, 'ts_map')
    analysis.ts_maps.append(ts_map_path)
    submit_fit(sargs, ts_map_path, sub_file=analysis.id+'.sub', ana_type='TS_Map', partition='xtralong')

    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --time_range {} {}'

    sargs = sargs.format(get_68_psf(ts_emin), analysis.fermi_data, analysis.ts_emin, analysis.ra, analysis.dec,
                         analysis.tsmjd1, analysis.tsmjd2)
    ts_map_short_path = os.path.join(analysis.bpath, 'ts_map_short')
    analysis.ts_maps.append(ts_map_short_path)
    submit_fit(sargs, ts_map_short_path, sub_file=analysis.id+'.sub', ana_type='TS_Map', partition='xtralong')

    #getsrcprob
    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
    sargs = sargs.format(get_68_psf(5000),analysis.fermi_data, 5000, analysis.ra, analysis.dec)
    srcprob_path = os.path.join(analysis.bpath, 'srcprob')
    analysis.srcprob_path = srcprob_path
    submit_fit(sargs, srcprob_path, srcs=analysis.srcs, sub_file=analysis.id+'.sub', ana_type='srcprob', partition='xtralong')

print('Submit_SEDs')

for src in analysis.srcs:
    if make_pdf:
        continue
    bpath_src = src_path(analysis.bpath, src.name)
    print(' \n \n {} is at a distance {:.1f} deg'.format(src.name, src.dist))
    if src.dist > args['max_dist']:
        continue
    if os.path.exists(bpath_src) and (not args['recovery']):
        print('Remove Path: {}'.format(bpath_src))
        shutil.rmtree(bpath_src)
    src.setup_folders(bpath_src)
    src.get_mw_data(this_path)
    
    src.make_sed(analysis.emin, analysis.fermi_data, name='', add_srcs=analysis.srcs, job_id = analysis.id)
    if analysis.emin < 1e3:
        src.make_sed(analysis.emin, analysis.fermi_data, name='1GeV', add_srcs=analysis.srcs,
                     job_id=analysis.id)

    src.make_fixed_binning_lightcurve(analysis.emin, analysis.fermi_data, analysis.mjd_range, mjd=analysis.mjd,
                                      dt_lc=args['dt_lc'], mode=args['mode'], name='', add_srcs=analysis.srcs,
                                      job_id = analysis.id)
    if analysis.emin < 1e3:
        src.make_fixed_binning_lightcurve(1e3, analysis.fermi_data, analysis.mjd_range, mjd=analysis.mjd,
                                          dt_lc=args['dt_lc'], mode=args['mode'], name='1GeV',
                                          add_srcs=analysis.srcs, job_id=analysis.id)

if os.path.exists(analysis_object_path):
    os.remove(analysis_object_path)
with open(analysis_object_path, "wb") as f:
    pickle.dump(analysis, f)
mins = 0
final_pdf = False
prev_len_jobs = -1

# Wait for jobs to finish....
while not final_pdf:
    if not make_pdf:
        time.sleep(60)
    mins += 1
    jobs = os.popen('squeue --user ga53lag').read()
    len_jobs = len([i for i in jobs.split('\n') if analysis.id in i])
    print len_jobs
    if len_jobs == 0 or make_pdf:
        final_pdf = True
    if (not mins % 60 == 1 or not len_jobs != prev_len_jobs) and not final_pdf:
        continue
    prev_len_jobs = len_jobs

## Make Plots
    if args['overwrite']:
        analysis.make_ts_maps()
        for src in analysis.srcs:
            if src.dist > args['max_dist']:
                continue
            src.make_lc_plot(analysis.mjd)
            if analysis.emin < 1e3: 
                src.make_sed_lightcurve(lcs=['default', '1GeV'])
            else: 
                src.make_sed_lightcurve(lcs=['default'])

    src_latex = ''
    for src in analysis.srcs:
        if src.dist > args['max_dist']:
            continue
        src_latex += src.source_summary(analysis.bpath, analysis.mjd, mode=args['mode'])
    analysis.make_pdf(src_latex, final_pdf = final_pdf)
    print_to_slack('Fit Results', analysis.pdf_out_path)
