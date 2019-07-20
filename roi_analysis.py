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
gcn_dict = None
rec = args['recovery']
make_pdf = args['make_pdf']
overwrite = args['overwrite']
print('Run with args')
print(args)

bpath = os.path.join(args['basepath'], args['event'])
if os.path.exists(os.path.join(bpath, 'analysis.pickle')):
    print('Folder {} already exist....exit. Only create PDF'.format(bpath))
    with open(os.path.join(bpath,'analysis.pickle'), "rb") as f:
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
    analysis.bpath = bpath
    analysis.emin = args['emin']
    analysis.event_name = ev_str
    if 'err90' in args.keys():
        analysis.err90 = args['err90']


this_path = os.path.dirname(os.path.abspath(__file__))

if not rec and not make_pdf:

    # Run VOU Tool
    analysis.ROI_analysis(this_path, args['radius'])
    if args['only_vou']:
        exit()

    # Download gamma-ray data
    analysis.get_fermi_data(days=args['dt'], mjd_range=args['mjd_range'])
    with open(os.path.join(analysis.bpath,'analysis.pickle'), "wb") as f:
        pickle.dump(analysis, f)

if not make_pdf:
    args['overwrite'] = True

# Start the gamma-ray analysis

# TS maps

    ts_emin = np.max([1000, analysis.emin])
    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
    sargs = sargs.format(get_68_psf(ts_emin), analysis.fermi_data, ts_emin, analysis.ra, analysis.dec)
    ts_map_path = os.path.join(analysis.bpath, 'ts_map')
    analysis.ts_maps.append(ts_map_path)
    submit_fit(sargs, ts_map_path, sub_file=analysis.id+'.sub', ana_type='TS_Map', partition='xtralong')

    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --time_range {} {}'
    if args['mode'] == 'end':
        tsmjd1 = analysis.mjd-200
        tsmjd2 = analysis.mjd

    else:
        tsmjd1 = analysis.mjd-100
        tsmjd2 = analysis.mjd+100
    sargs = sargs.format(get_68_psf(ts_emin), analysis.fermi_data, ts_emin, analysis.ra, analysis.dec,
                         tsmjd1, tsmjd2)
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
    bpath_src = src_path(analysis.bpath, src.name)
    print('{} is at a distance {:.1f} deg'.format(src.name, src.dist))
    if src.dist > args['max_dist']:
        continue
    if os.path.exists(bpath_src) and (not args['recovery']) and (not args['make_pdf']):
        print('Remove Path: {}'.format(bpath_src))
        shutil.rmtree(bpath_src)
    src.setup_folders(bpath_src)
    if make_pdf:
        continue
    src.get_mw_data(this_path)
    if make_pdf:
        continue
    
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

if os.path.exists(os.path.join(analysis.bpath,'analysis.pickle')):
    os.remove(os.path.join(analysis.bpath,'analysis.pickle'))
with open(os.path.join(analysis.bpath,'analysis.pickle'), "wb") as f:
    pickle.dump(analysis, f)
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

    if args['overwrite']:
        for tsm in analysis.ts_maps:
            try:
                plot.make_ts_plot(tsm, analysis.srcs,
                                  os.path.join(analysis.vou_out, 'find_out_temp.txt'),
                                  plt_mode='tsmap', error90 = analysis.err90)
            except Exception as inst:
                warnings.warn("Couldn't create ts map...")
                print(inst)

            try:
                plot.make_ts_plot(tsm, analysis.srcs,
                                  os.path.join(analysis.vou_out, 'find_out_temp.txt'),
                                  plt_mode='residmap',  error90=analysis.err90)
            except Exception as inst:
                warnings.warn("Couldn't create residual map...")
                print(inst)

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
    with open(os.path.join(analysis.bpath, 'vou_blazar/full_output'), 'r') as f:
        full_out = f.read().encode('utf8').\
            replace('\x1b', '').replace('nu_p', 'nu peak')
        full_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', full_out)
        full_out = re.sub('(Match[^\\\\]*)', r'\\textbf{\1}', full_out)
        full_out = re.sub('(Dist[^\\\\]*)', r'\\textbf{\1}', full_out)
        full_out = [i.strip() + ' \n' for i in full_out.split('\n')]
        full_out = [i if i != '\\\\ \n' else '\\\\\\\ \n' for i in full_out]
        full_out = ' '.join(full_out)
    with open(os.path.join(analysis.bpath, 'vou_blazar/short_output'), 'r') as f:
        short_out = f.read().encode('utf8').replace('\x1b', '')
        short_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', short_out)
    with open('./template.tex') as f:
        template = f.read()

    c = SkyCoord(analysis.ra, analysis.dec, frame='icrs', unit="deg")
    gal = c.galactic
    if not final_pdf:
        prelim = '{\color{red} Warning: The following results are all preliminary and will be continously updated whenever calculations are finished.}'
    else:
        prelim = ' '
    t = Time.now()
    t_now_str = t.iso
    out = template.format(prelim=prelim,
                          date=args['mjd'],
                          ra=analysis.ra,
                          dec=analysis.dec,
                          emin=1.*analysis.emin/1000.,
                          l=gal.l.deg,
                          b=gal.b.deg,
                          cat_srcs=short_out,
                          rx_map=os.path.join(analysis.vou_out, 'RX_map.eps'),
                          vou_pic=os.path.join(analysis.vou_out, 'candidates.eps'),
                          ts_map=os.path.join(analysis.bpath, 'ts_map/tsmap.png'),
                          ts_map_short=os.path.join(analysis.bpath, 'ts_map_short/tsmap.png'),
                          vou_output=full_out,
                          event=analysis.event_name,
                          tsmjd1=tsmjd1,
                          tsmjd2=tsmjd2,
                          src_latex=src_latex,
                          mjd1=analysis.mjd_range[0],
                          mjd2=analysis.mjd_range[1],
                          energy=analysis.emin/1000.,
                          tsemin = ts_emin/1e3,
                          gcnurl = analysis.gcn,
                          createdon=t_now_str)
    latex_path = os.path.join(analysis.bpath, analysis.event_name + '.tex')
    if os.path.exists(latex_path):
        os.remove(latex_path)
    with open(latex_path, 'w+') as f:
        f.write(out)
    if not os.path.exists(os.path.join(analysis.bpath, 'sample.bib')):
        shutil.copyfile('sample.bib', os.path.join(analysis.bpath, 'sample.bib'))
    os.chdir(analysis.bpath)
    cmds  = [
        ['pdflatex', '-interaction', 'nonstopmode', analysis.event_name + '.tex'],
        ['bibtex', analysis.event_name + '.aux'],
        ['pdflatex', '-interaction', 'nonstopmode',  analysis.event_name + '.tex'],
        ['pdflatex', '-interaction', 'nonstopmode',  analysis.event_name + '.tex'],
    ]
    for c in cmds:
        subprocess.call(c)
    os.chdir(this_path)
    print_to_slack('Fit Results', os.path.join(analysis.bpath, ev_str + '.pdf'))
