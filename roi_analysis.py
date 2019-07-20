# coding: utf-8

import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('/scratch9/tglauch/realtime_service/python/lib/python2.7/site-packages')
sys.path.append('/home/ga53lag/Software/python_scripts/')
from get_fermi_data import get_data
import subprocess
import argparse
import os
from roi_functions import get_lc_time3fgl, get_lc_time4fgl, make_gif, submit_fit, \
                          path_settings, vou_path, partition_t, get_68_psf 
import shutil
import datetime
from slack_lib import print_to_slack
import numpy as np
from myfunctions import MET_to_MJD, dict_to_nparray, ts_to_pval, pval_to_sigma
import plot
import warnings
import re
import time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
from collections import OrderedDict
from gcn import read_gcn
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
        required=False)
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

if os.path.exists(os.path.join(bpath, 'run_info.npy')) and not rec and not args['only_vou']:
    print('Folder already exist....exit. Only create PDF')
    make_pdf=True

## read GCN
if args['gcn'] is not None:
    analysis = Analysis()
    analysis.read_gcn(args['gcn'])
else:
    if args['mjd'] is not None:
        t = Time(args['mjd'], format='mjd')
    else:
        t = Time.now()
        args['mjd'] = t.mjd
    analysis = Analysis(args['mjd'], args['ra'], args['dec'])

if args['event'] is not None:
    ev_str = args['event']
else:
    t = Time(analysis.mjd, format='mjd')
    ev_str = 'IC{}{:02d}{:02d}'.format(str(t.datetime.year)[-2:],
                                       t.datetime.month, t.datetime.day)
bpath = os.path.join(args['basepath'], ev_str)
analysis.bpath = bpath
analysis.emin = args['emin']
analysis.event_name = ev_str
analysis.fermi_data = os.path.join(bpath, 'fermi_data')
analysis.vou_out = os.path.join(bpath, 'vou_blazar')
if 'err90' in args.keys():
    analysis.err90 = args['err90']


this_path = os.path.dirname(os.path.abspath(__file__))

if not rec and not make_pdf:
    if not os.path.exists(analysis.vou_out):
        os.makedirs(analysis.vou_out)
    os.chdir(analysis.vou_out)
    cmd = [vou_path,
           str(analysis.ra), str(analysis.dec), args['radius'], str(30), str(90)]
    for i in range(2):
        subprocess.call(cmd)

    # Setup Variables
    out_str = analysis.get_sources()
    headline = '*Result for {} *\n'.format(ev_str)
    sum_text = headline + out_str
    print_to_slack(sum_text)

    # Convert VOU Blazar Output
    rx_ps = os.path.join(analysis.vou_out, 'RX_map.ps')
    os.system('ps2eps -B ' + rx_ps)
    cand_ps = os.path.join(analysis.vou_out, 'candidates.ps')
    os.system('ps2eps -B ' + cand_ps )

    #Create VOU Source Summary
    with open('./short_output', 'w+') as ofile:
        ofile.write(out_str.replace('\n' ,'\\\ \n'))
    with open('./phase1', 'r') as ifile:
        lines = ifile.read().split('Gamma-ray Counterparts')[0]
        lines = re.sub('\\[..?;.?m', ' ', lines)
        lines = lines.replace('[0m', ' ')
        print_to_slack('', pic=os.path.join(analysis.vou_out, 'candidates.eps'))
    with open('./full_output', 'w+') as ofile:
        ofile.write(lines.replace('\n' ,'\\\ \n'))
    os.chdir(this_path)
    if args['only_vou']:
        exit()

    # download gamma-ray data
    kwargs = {'emin': analysis.emin, 'out_dir' : analysis.fermi_data}
    if 'dt' in args.keys():
        kwargs['days'] = args['dt']
    if args['mjd_range'] is None:
        kwargs['mjd'] = args['mjd']
        kwargs['mode'] = args['mode']
    else:
        kwargs['mjd'] = args['mjd_range']
    print kwargs
    MET = get_data(args['ra'], args['dec'], **kwargs)
    MJD = [MET_to_MJD(float(i)) for i in MET]
    analysis.mjd_range = MJD
    with open(os.path.join(bpath,'analysis.pickle'), "wb") as f:
        pickle.dump(analysis, f)
else:
    run_info = np.load(os.path.join(bpath,'run_info.npy'), allow_pickle=True)[()]
    with open(os.path.join(bpath,'analysis.pickle'), "rb") as f:
        analysis = pickle.load(f)
    args['recovery'] = rec
    args['make_pdf'] = make_pdf
    args['overwrite'] = overwrite
    if 'gcn' in args.keys():
        print('Update GCN Info from {}'.format(args['gcn']))
        read_gcn(args['gcn'], bpath)
        analysis.ra = gcn_dict['SRC_RA']
        analysis.dec = gcn_dict['SRC_DEC']
        analysis.err90 = gcn_dict['SRC_ERROR']
        analysis.err50 = gcn_dict['SRC_ERROR50']

if not make_pdf:
    args['overwrite'] = True
# start the gamma-ray analysis
#TS maps
    ts_emin = np.max([1000, analysis.emin])
    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
    sargs = sargs.format(get_68_psf(ts_emin), analysis.fermi_data, ts_emin, analysis.ra, analysis.dec)
    ts_map_path = os.path.join(bpath, 'ts_map')
    submit_fit(sargs, ts_map_path, sub_file=ev_str+'.sub', ana_type='TS_Map', partition='xtralong')

    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --time_range {} {}'
    if args['mode'] == 'end':
        tsmjd1 = analysis.mjd-200
        tsmjd2 = analysis.mjd

    else:
        tsmjd1 = analysis.mjd-100
        tsmjd2 = analysis.mjd+100
    sargs = sargs.format(get_68_psf(ts_emin), analysis.fermi_data, ts_emin, analysis.ra, analysis.dec,
                         tsmjd1, tsmjd2)
    ts_map_short_path = os.path.join(bpath, 'ts_map_short')
    submit_fit(sargs, ts_map_short_path, sub_file=ev_str+'.sub', ana_type='TS_Map', partition='xtralong')

    #getsrcprob
    sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
    sargs = sargs.format(get_68_psf(5000),analysis.fermi_data, 5000, analysis.ra, analysis.dec)
    srcprob_path = os.path.join(bpath, 'srcprob')
    submit_fit(sargs, srcprob_path, srcs=analysis.srcs, sub_file=ev_str+'.sub', ana_type='srcprob', partition='xtralong')

print('Submit_SEDs')

for src in analysis.srcs:
    bpath_src = src_path(analysis.bpath, src.name)
    print('{} is at a distance {}'.format(src.name, src.dist))
    if src.dist > args['max_dist']:
        print('Source exceeds distance of {} deg. No job will be submitted'.format(args['max_dist']))
        continue
    src.setup_folders(bpath_src)
    if os.path.exists(bpath_src) and (not args['recovery']) and (not args['make_pdf']):
        print('Remove Path: {}'.format(bpath_src))
        shutil.rmtree(bpath_src)
    if make_pdf:
        continue
    src.get_mw_data(this_path)
    if make_pdf:
        continue
    
    src.make_sed(analysis.emin, analysis.fermi_data, name='', add_srcs=analysis.srcs)
    if analysis.emin < 1e3:
        src.make_sed(analysis.emin, analysis.fermi_data, name='1GeV', add_srcs=analysis.srcs)

    src.make_fixed_binning_lightcurve(analysis.emin, analysis.fermi_data, analysis.mjd_range, mjd=analysis.mjd,
                                      dt_lc=args['dt_lc'], mode=args['mode'], name='', add_srcs=analysis.srcs)
    if analysis.emin < 1e3:
        src.make_fixed_binning_lightcurve(1e3, analysis.fermi_data, analysis.mjd_range, mjd=analysis.mjd,
                                          dt_lc=args['dt_lc'], mode=args['mode'], name='1GeV',
                                          add_srcs=analysis.srcs)
mins = 0
final_pdf = False
prev_len_jobs = -1
while not final_pdf:
    if not make_pdf:
        time.sleep(60)
    mins += 1
    jobs = os.popen('squeue --user ga53lag').read()
    len_jobs = len([i for i in jobs.split('\n') if ev_str in i])
    print len_jobs
    if len_jobs == 0 or make_pdf:
        final_pdf = True
    if (not mins % 60 == 1 or not len_jobs != prev_len_jobs) and not final_pdf:
        continue
    prev_len_jobs = len_jobs

    if args['overwrite']:
        try:
            if final_pdf:
                yaxis = False
            else:
                yaxis = True
            plot.make_ts_plot(ts_map_path, analysis.srcs,
                              os.path.join(analysis.vou_out, 'find_out_temp.txt'),
                              plt_mode='tsmap', yaxis=yaxis, error90 = analysis.err90)
        except Exception as inst:
            warnings.warn("Couldn't create ts map...")
            print(inst)

        try:
            plot.make_ts_plot(ts_map_short_path, analysis.srcs,
                              os.path.join(analysis.vou_out, 'find_out_temp.txt'),
                              plt_mode='tsmap', legend=False, error90= analysis.err90)
        except Exception as inst:
            warnings.warn("Could't create residual map...")
            print(inst)

        try:
            plot.make_ts_plot(ts_map_path, analysis.srcs,
                              os.path.join(analysis.vou_out, 'find_out_temp.txt'),
                              plt_mode='residmap', legend=False, error90=analysis.err90)
        except Exception as inst:
            warnings.warn("Couldn't create residual map...")
            print(inst)

        for src in analysis.srcs:
            src.make_lc_plot(analysis.mjd)
            if analysis.emin < 1e3: 
                src.make_sed_lightcurve(lcs=['default', '1GeV'])
            else: 
                src.make_sed_lightcurve(lcs=['default'])

    src_latex = ''
    for src in analysis.srcs:
        if src.dist > args['max_dist']:
            print('Source exceeds distance of {} deg. No summary information will be added'.format(args['max_dist']))
            continue
        src_latex += src.source_summary(bpath, analysis.mjd, mode=args['mode'])
    with open(os.path.join(bpath, 'vou_blazar/full_output'), 'r') as f:
        full_out = f.read().encode('utf8').\
            replace('\x1b', '').replace('nu_p', 'nu peak')
        full_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', full_out)
        full_out = re.sub('(Match[^\\\\]*)', r'\\textbf{\1}', full_out)
        full_out = re.sub('(Dist[^\\\\]*)', r'\\textbf{\1}', full_out)
        full_out = [i.strip() + ' \n' for i in full_out.split('\n')]
        full_out = [i if i != '\\\\ \n' else '\\\\\\\ \n' for i in full_out]
        full_out = ' '.join(full_out)
    with open(os.path.join(bpath, 'vou_blazar/short_output'), 'r') as f:
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
    gcnurl = 'not yet implemented'
    if 'gcn' in args.keys():
        gcnurl = args['keys']
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
                          ts_map=os.path.join(bpath, 'ts_map/tsmap.png'),
                          ts_map_short=os.path.join(bpath, 'ts_map_short/tsmap.png'),
                          vou_output=full_out,
                          event=ev_str,
                          tsmjd1=tsmjd1,
                          tsmjd2=tsmjd2,
                          src_latex=src_latex,
                          mjd1=analysis.mjd_range[0],
                          mjd2=analysis.mjd_range[1],
                          energy=analysis.emin/1000.,
                          tsemin = ts_emin/1e3,
                          gcnurl = gcnurl,
                          createdon=t_now_str)
    latex_path = os.path.join(bpath,ev_str + '.tex')
    if os.path.exists(latex_path):
        os.remove(latex_path)
    with open(latex_path, 'w+') as f:
        f.write(out)
    if not os.path.exists(os.path.join(bpath, 'sample.bib')):
        shutil.copyfile('sample.bib', os.path.join(bpath, 'sample.bib'))
    os.chdir(bpath)
    cmds  = [
        ['pdflatex', '-interaction', 'nonstopmode', ev_str + '.tex'],
        ['bibtex', ev_str + '.aux'],
        ['pdflatex', '-interaction', 'nonstopmode',  ev_str + '.tex'],
        ['pdflatex', '-interaction', 'nonstopmode',  ev_str + '.tex'],
    ]
    for c in cmds:
        subprocess.call(c)
    os.chdir(this_path)
    print_to_slack('Fit Results', os.path.join(bpath, ev_str + '.pdf'))
