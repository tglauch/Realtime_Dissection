# coding: utf-8


# run for example as python roi_analysis.py --mjd 55497.30 --ra 88.95 --dec 0.50 --adaptive_scaling --max_dist 0.8 --err90 0.48 -0.53 0.25 -0.21
import logging
#logging.getLogger().setLevel(logging.INFO)
logging.info('Imports')
import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('/home/ga53lag/Software/python_scripts/')
sys.path.append('./lib')
import argparse
import os
import datetime
import numpy as np
import re
import time
#import plot
from astropy.time import Time
from analysis import Analysis
import pickle
from slack_lib import print_to_slack
from add_classes import Ellipse


def parseArguments():
    """Parse the command line arguments
    Returns:
    args : Dictionary containing the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ra", help="right ascension",type=float)
    parser.add_argument(
        "--dec", help="declination", type=float)
    parser.add_argument(
        "--mjd", help="The MJD time of the event", type=float)
    parser.add_argument(
        "--mjd_range", help="The mjd range for the analysis", type=float, nargs='+')
    parser.add_argument(
        "--dt_lc", help="Time lenght of bins in the light curve", type=float, default=200)
    parser.add_argument(
        "--dt", help="Length of time window to analyze", type=float)
    parser.add_argument(
        "--radius", help="Radius of the region to analyze", type=str, default="180")
    parser.add_argument(
        "--emin", help="Lower energy bound for SED", type=float, default=100)
    parser.add_argument(
        "--event", help="Name of the Event", default='None')
    parser.add_argument(
        "--vou", help="Only run VOU Blazar", action="store_true", default = False)
    parser.add_argument(
        "--gcn", help="Pass a link to a GCN notice", type=str)
    parser.add_argument(
        "--mode", help="end if given mjd is at the end of the time window, mid if in the middle'",
        type=str, default='end')
    parser.add_argument(
        "--lat_analysis", help="Run the LAT analysis?",
        action="store_true", default = False)
    parser.add_argument(
        "--make_pdf", help="Only create the output PDF form a previous run ",
        action="store_true", default = False)
    parser.add_argument(
        "--overwrite", help="Overwrite already existing Plots",
        action="store_true", default = False)
    parser.add_argument(
        "--basepath", help = 'Basepath for the Output',
        type=str, default = '/scratch9/tglauch/realtime_service/output/')
    parser.add_argument(
        "--max_dist", help = 'Radius of sources to be included', type=float, default=2.5)
    parser.add_argument(
        "--err90", help= 'The 90 percent error ellipse. Format:  ra1 ra2 dec1 dec2', nargs='+')
    parser.add_argument(
        "--adaptive_scaling", help="Scale VOU and TS map region adaptively",
        action="store_true", default = False) 
    args = parser.parse_args()
    return args.__dict__


args = parseArguments()
make_pdf = args['make_pdf']
print('Run with args \n {}'.format(args))

# Setup Analysis Class
bpath = os.path.join(args['basepath'], args['event'])
analysis_object_path = os.path.join(bpath,'analysis.pickle')
if os.path.exists(analysis_object_path):
    logging.info('Read analysis object from {}'.format(bpath))
    with open(analysis_object_path, "rb") as f:
        analysis = pickle.load(f)
    analysis.update_gcn()
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
    analysis.emin = float(args['emin'])
    analysis.event_name = ev_str
    analysis.radius = float(args['radius'])
    analysis.max_dist = float(args['max_dist'])
    args['vou'] = True
    args['lat_analysis'] = True
analysis.this_path = os.path.dirname(os.path.abspath(__file__))
if args['err90'] is not None:
    if len(args['err90']) == 1:
        analysis.err90 = float(args['err90'])
    else:
        analysis.err90 = Ellipse(analysis.ra, analysis.dec, args['err90'])
if args['adaptive_scaling']:
        analysis.adaptive_radius()
if args['vou']:
    # Run VOU Tool
    logging.info('Start the VOU analysis pipeline')
    analysis.ROI_analysis()
    if os.path.exists(analysis_object_path):
        os.remove(analysis_object_path)
    with open(analysis_object_path, "wb") as f:
        pickle.dump(analysis, f)

# Start the gamma-ray analysis
if args['lat_analysis']:
    logging.info('Start the Fermi LAT Analysis pipeline')
    args['overwrite'] = True
    # Download gamma-ray data
    analysis.get_fermi_data(days=args['dt'], mjd_range=args['mjd_range'])
    with open(analysis_object_path, "wb") as f:
        pickle.dump(analysis, f)

    if analysis.mode == 'end':
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
    ts_map_path = os.path.join(analysis.bpath, 'ts_map') 
    analysis.make_ts_map(ts_map_path) 
    ts_map_short_path = os.path.join(analysis.bpath, 'ts_map_short') 
    analysis.make_ts_map(ts_map_short_path, trange=[analysis.tsmjd1, analysis.tsmjd2])
    ts_map_short_path = os.path.join(analysis.bpath, 'ts_map_very_short')
    analysis.make_ts_map(ts_map_short_path, trange=[analysis.mjd-1, analysis.mjd+1])
    #Calculate Source Probability
    srcprob_path = os.path.join(analysis.bpath, 'srcprob')
    analysis.calc_src_probs(srcprob_path, emin=5000)

    # Run SEDs
    print('Submit_SEDs')
    for src in analysis.srcs:
        if make_pdf:
            continue
        bpath_src = os.path.join(bpath, src.name.replace(' ', '_'))
        print(' \n \n {} is at a distance {:.1f} deg'.format(src.name, src.dist))
        if src.dist > analysis.max_dist:
            continue
        src.setup_folders(bpath_src)
        src.get_mw_data(analysis.this_path)
        add_srcs = [src] #analysis.srcs 
        src.make_sed(analysis.emin, analysis.fermi_data, name='', add_srcs=add_srcs, job_id = analysis.id)
        if analysis.emin < 1e3:
            src.make_sed(analysis.emin, analysis.fermi_data, name='1GeV', add_srcs=add_srcs,
                         job_id=analysis.id)

        src.make_fixed_binning_lightcurve(analysis.emin, analysis.fermi_data, analysis.mjd_range, mjd=analysis.mjd,
                                          dt_lc=args['dt_lc'], mode=analysis.mode, name='', add_srcs=add_srcs,
                                          job_id = analysis.id)
        if analysis.emin < 1e3:
            src.make_fixed_binning_lightcurve(1e3, analysis.fermi_data, analysis.mjd_range, mjd=analysis.mjd,
                                              dt_lc=args['dt_lc'], mode=analysis.mode, name='1GeV',
                                              add_srcs=add_srcs, job_id=analysis.id)

    if os.path.exists(analysis_object_path):
        os.remove(analysis_object_path)
    with open(analysis_object_path, "wb") as f:
        pickle.dump(analysis, f)

# Wait for jobs to finish....
mins = 0
final_pdf = False
prev_len_jobs = -1
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

    # Make Plots
    logging.info('Creating Plots...')
    if args['overwrite']:
        analysis.make_ts_map_plots()
        for src in analysis.srcs:
            if src.dist > analysis.max_dist:
                continue
            src.make_lc_plot(analysis.mjd)
            if analysis.emin < 1e3: 
                src.make_sed_lightcurve(lcs=['default', '1GeV'])
            else: 
                src.make_sed_lightcurve(lcs=['default'])
    logging.info('Make PDF..')
    src_latex = ''
    for src in analysis.srcs:
        if src.dist > analysis.max_dist:
            continue
        src_latex += src.source_summary(analysis.bpath, analysis.mjd, mode=analysis.mode)
    analysis.create_html()
    analysis.make_pdf(src_latex, final_pdf = final_pdf)
    print_to_slack('Fit Results', analysis.pdf_out_path)
