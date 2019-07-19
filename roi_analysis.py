# coding: utf-8

import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('/scratch9/tglauch/realtime_service/python/lib/python2.7/site-packages')
sys.path.append('/home/ga53lag/Software/python_scripts/')
from get_fermi_data import get_data
import subprocess
import argparse
import os
from roi_functions import get_sources, get_lc_time3fgl, get_lc_time4fgl, make_gif, submit_fit, generate_src_xml, \
                          path_settings, vou_path, partition_t 
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
from scipy.interpolate import interp1d
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


def get_68_psf(E):
    x = np.genfromtxt('./lat_68_psf.txt', delimiter = ',')
    return interp1d(x[:,0], x[:,1])(E)


def source_summary(bpath, src, mjd, mode='mid'):
    bpath_src = src_path(bpath, src['name'])
    lc_base = os.path.join(bpath_src, path_settings['lc'])
    lc_path = os.path.join(lc_base, 'lightcurve.pdf')
    folders = [fold for fold in os.listdir(lc_base) if os.path.isdir(os.path.join(lc_base, fold))]
    if len(folders) == 0:
        return ''
    if mode == 'end':
        print sorted(folders)
        sed_path = os.path.join(lc_base,sorted(folders)[-1])
    else:
        for fold in folders:
            if (mjd <= float(fold.split('_')[1])) and (mjd > float(fold.split('_')[0])):
                sed_path = os.path.join(lc_base,fold)
                break
    sed_full_res = np.load(os.path.join(bpath_src, path_settings['sed'], 'llh.npy'), allow_pickle=True)[()]
    ts = sed_full_res['sources'][src['name']]['ts']
    sigma = np.max([0, pval_to_sigma(ts_to_pval(ts, 1))])
    l_str ='\subsection{{{srcinfo}}}'
    if os.path.exists(os.path.join(sed_path, 'llh.npy')):
        fit_res = np.load(os.path.join(sed_path, 'llh.npy'), allow_pickle=True)[()]
        ra = fit_res['sources'][src['name']]['RAJ2000']
        dec = fit_res['sources'][src['name']]['DEJ2000']
        t_str = '{src_name} $|$ \\small\\textnormal{{{src_info}}}'
        t_str_info = 'ra = {ra:.2f}$^\circ$, dec = {dec:.2f}$^\circ$, $\Sigma$= {sigma:.1f} $\sigma$, $\Delta\psi$ = {dist:.2f}$^\circ$'
        t_str_info = t_str_info.format(sigma=sigma, ra=ra, dec=dec, dist=src['dist'])
        t_str = t_str.format(src_name = src['name'], src_info=t_str_info)
        l_str = l_str.format(srcinfo=t_str) # srcname=src['name'], srcinfo=t_str)
    else:
        t_str = '{src_name} $|$ \\small\\textnormal{{{src_info}}}'
        t_str_info = 'ra = {ra:.2f}$^\circ$, dec = {dec:.2f}$^\circ$, $\Delta\psi$ = {dist:.2f}$^\circ$'
        t_str_info = t_str_info.format(ra=src['ra'], dec=src['dec'], dist=src['dist'])
        t_str = t_str.format(src_name = src['name'], src_info=t_str_info)
        l_str = l_str.format(srcinfo=t_str)
    print src
    l_str += '{} \\'.format(src['alt_name'])
    sed_path = os.path.join(sed_path, 'sed.pdf')
    try:
        srcprob = fits.open(os.path.join(bpath,'srcprob/ft1_srcprob_00.fits'))
        energy = srcprob[1].data['ENERGY']
        mjd = MET_to_MJD(srcprob[1].data['TIME'])
        src_prob = srcprob[1].data[src['name']]
        prob_mask = src_prob > 0.90
        ind = np.argsort(energy[prob_mask])[::-1]
        with open('./tab_template.tex', 'r') as f:
            tab_str = f.read()
        prob_str = ''
        for i in range(np.min([5, len(mjd[prob_mask][ind])])):
            prob_str += '{:.2f} & {:.2f} & {:.2f} \\\\ \n'.format(mjd[prob_mask][ind[i]],
                                                                  src_prob[prob_mask][ind[i]]*100,
                                                                  energy[prob_mask][ind[i]]/1e3)
        tab_str = tab_str.format(events=prob_str)
        l_str += tab_str
    except Exception as inst:
        warnings.warn("Could not find source probabilities")
        print(inst)
    print l_str
    with open('source_summary.tex', 'r') as infile:
        fig_str = infile.read()       
    if os.path.exists(sed_path):
        cap_str = 'SED for {}. See the description in section \\ref{{sec:sed}} for more details.'
        l_str += fig_str.format(width = 0.7, path = sed_path, caption=cap_str.format(src['name']))
    if os.path.exists(lc_path):
        l_str += fig_str.format(width = 0.7, path = lc_path, caption='Light curve for {}'.format(src['name']))
    gev_lc = os.path.join(lc_base+'_1GeV', 'lightcurve.pdf')
    if os.path.exists(gev_lc):
        l_str += fig_str.format(width = 0.8, path = gev_lc, caption='1GeV light curve for {}'.format(src['name']))
    if ('4FGL' in src['name']) | ('3FGL' in src['name']):
        cp_candidate_path = os.path.join(bpath_src, 'vou_counterpart', src['name'].replace(' ', '_').replace('.','_') + '.eps') 
        if os.path.exists(cp_candidate_path):
            l_str += fig_str.format(width = 0.7, path = cp_candidate_path, caption='Possible couterparts for {}'.format(src['name']))
    l_str += '\\clearpage \n'
    return l_str


def make_fixed_binning_lc(src, sub_args,  dt_lc, mjd_range, opath,
                          sub_file='lc.sub', mjd_mid=None, mode='end', **args):
    if args['mode'] == 'end':
        time_windows = [[k - dt_lc, k] for k in
                         np.abs(np.arange(-mjd_range[1], -mjd_range[0], dt_lc))]
    elif args['mode'] == 'mid':
        time_windows = [[k - dt_lc, k] for k in
                         np.abs(np.arange(mjd_mid + 3 * dt_lc / 2, mjd_range[1], dt_lc))]
        time_windows2 = [[k - dt_lc, k] for k in
                         np.abs(np.arange(-mjd_mid - dt_lc / 2, -mjd_range[0], dt_lc))]
        time_windows.extend(time_windows2)

    if (time_windows[-1][1] - mjd_range[0]) < dt_lc:
        del time_windows[-1]
    time_windows.append('')

    for t_window in time_windows:
        if t_window == '':
            partition='long'
        else:
            add_str = '{:.1f}_{:.1f}'.format(t_window[0], t_window[1])
            partition='kta'
        submit_fit(sub_args, opath, src, sub_file=sub_file, trange=t_window, partition=partition, **args)
    return


def src_path(bpath, src):
    return os.path.join(bpath, src.replace(' ', '_'))


args = parseArguments()
gcn_dict = None
rec = args['recovery']
make_pdf = args['make_pdf']
overwrite = args['overwrite']
print('Run with args')
print(args)
if args['gcn'] is not None:
    gcn_dict = read_gcn(args['gcn'])
    args['ra'] = gcn_dict['SRC_RA']
    args['dec'] = gcn_dict['SRC_DEC']
    args['err90'] = gcn_dict['SRC_ERROR']
    args['err50'] = gcn_dict['SRC_ERROR50']
    args['mjd'] = gcn_dict['MJD']
    t = Time.now()
    if np.abs(args['mjd'] - t.mjd)<2:
        args['mode']= 'end'
    else:
        args['mode'] = 'mid'
if args['mjd'] is not None:
    t = Time(args['mjd'], format='mjd')
else:
    t = Time.now()
    args['mjd'] = t.mjd
if args['event'] is not None:
    ev_str = args['event']
else:
    ev_str = 'IC{}{:02d}{:02d}'.format(str(t.datetime.year)[-2:],
                                       t.datetime.month, t.datetime.day)
bpath = os.path.join(args['basepath'], ev_str)
analysis = Analysis(ev_str, args['ra'], args['dec'])
analysis.bpath = bpath
if os.path.exists(os.path.join(bpath, 'run_info.npy')) and not rec and not args['only_vou']:
    print('Folder already exist....exit. Only create PDF')
    make_pdf=True
this_path = os.path.dirname(os.path.abspath(__file__))
fermi_data = os.path.join(bpath, 'fermi_data')
vou_out = os.path.join(bpath, 'vou_blazar')

if not rec and not make_pdf:
    if not os.path.exists(vou_out):
        os.makedirs(vou_out)
    if os.path.exists(os.path.join(bpath,'run_info.npy')):
        run_info = np.load(os.path.join(bpath,'run_info.npy'), allow_pickle=True)[()]
        args['ra'] = run_info['args']['ra']
        args['dec'] = run_info['args']['dec']
    os.chdir(vou_out)
    cmd = [vou_path,
           str(args['ra']), str(args['dec']), args['radius'], str(30), str(90)]
    for i in range(2):
        subprocess.call(cmd)

    # Setup Variables
    out_str = analysis.get_sources()
    headline = '*Result for {} *\n'.format(ev_str)
    sum_text = headline + out_str
    print_to_slack(sum_text)

    # Convert VOU Blazar Output
    rx_ps = os.path.join(vou_out, 'RX_map.ps')
    os.system('ps2eps -B ' + rx_ps)
    cand_ps = os.path.join(vou_out, 'candidates.ps')
    os.system('ps2eps -B ' + cand_ps )

    #Create VOU Source Summary
    with open('./short_output', 'w+') as ofile:
        ofile.write(out_str.replace('\n' ,'\\\ \n'))
    with open('./phase1', 'r') as ifile:
        lines = ifile.read().split('Gamma-ray Counterparts')[0]
        lines = re.sub('\\[..?;.?m', ' ', lines)
        lines = lines.replace('[0m', ' ')
        print_to_slack('', pic=os.path.join(vou_out, 'candidates.eps'))
    with open('./full_output', 'w+') as ofile:
        ofile.write(lines.replace('\n' ,'\\\ \n'))
    os.chdir(this_path)
    if args['only_vou']:
        exit()

    # download gamma-ray data
    kwargs = {'emin': args['emin'], 'out_dir' : fermi_data}
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
    #run_info = {'MJD' : MJD,
    #            'args' : args,
    #            'src_dict' : src_dict}
    #np.save(os.path.join(bpath,'run_info.npy'), run_info)
    with open(os.path.join(bpath,'analysis.pickle'), "wb") as f:
        pickle.dump(analysis, f)
else:
    run_info = np.load(os.path.join(bpath,'run_info.npy'), allow_pickle=True)[()]
    MJD = run_info['MJD']
    args = run_info['args']
    args['recovery'] = rec
    args['make_pdf'] = make_pdf
    args['overwrite'] = overwrite
    if 'gcn' in args.keys():
        print('Update GCN Info from {}'.format(args['gcn']))
        read_gcn(args['gcn'], bpath)
        args['ra'] = gcn_dict['SRC_RA']
        args['dec'] = gcn_dict['SRC_DEC']
        args['err90'] = gcn_dict['SRC_ERROR']
        args['err50'] = gcn_dict['SRC_ERROR50']
    #src_dict = run_info['src_dict']
#print src_dict
if not make_pdf:
    args['overwrite'] = True

# start the gamma-ray analysis

#TS maps
ts_emin = np.max([1000, args['emin']])
sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
sargs = sargs.format(get_68_psf(ts_emin), fermi_data, ts_emin, args['ra'], args['dec'])
ts_map_path = os.path.join(bpath, 'ts_map')
submit_fit(sargs, ts_map_path, sub_file=ev_str+'.sub', ana_type='TS_Map', partition='xtralong', **args)

sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --time_range {} {}'
if args['mode'] == 'end':
    tsmjd1 = args['mjd']-200
    tsmjd2 = args['mjd']

else:
    tsmjd1 = args['mjd']-100
    tsmjd2 = args['mjd']+100
sargs = sargs.format(get_68_psf(ts_emin), fermi_data, ts_emin, args['ra'], args['dec'],
                     tsmjd1, tsmjd2)
ts_map_short_path = os.path.join(bpath, 'ts_map_short')
submit_fit(sargs, ts_map_short_path, sub_file=ev_str+'.sub', ana_type='TS_Map', partition='xtralong', **args)

#getsrcprob
sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {}'
sargs = sargs.format(get_68_psf(5000),fermi_data, 5000, args['ra'], args['dec'])
srcprob_path = os.path.join(bpath, 'srcprob')
submit_fit(sargs, srcprob_path, srcs=analysis.srcs, sub_file=ev_str+'.sub', ana_type='srcprob', partition='xtralong', **args)

print('Submit_SEDs')
if 'max_dist' not in args.keys():
    args['max_dist'] = 2.

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
    src.get_mw_data()
    if make_pdf:
        continue
    sargs = '--target_src {} --free_radius {} --data_path {} --use_4FGL --emin {} '
    sub_args = sargs.format(src['name'].replace(' ', '_'), get_68_psf(args['emin']), fermi_data, args['emin'])
    tsub_args = sargs.format(src['name'].replace(' ', '_'), get_68_psf(1e3),  fermi_data, 1e3)
    if '3FGL' in src['name']:
        dt_lc = get_lc_time3fgl(src['name'], emin=args['emin'] ) 
    elif '4FGL' in src['name']:
        dt_lc = get_lc_time4fgl(src['name'], emin=args['emin'] ) 
    else:
        dt_lc = args['dt_lc']

    if args['mode'] == 'end':
        time_windows = [[k - dt_lc, k] for k in
                         np.abs(np.arange(-MJD[1], -MJD[0], dt_lc))]
    elif args['mode'] == 'mid':
        time_windows = [[k - dt_lc, k] for k in
                         np.abs(np.arange(args['mjd']+3*dt_lc/2, MJD[1], dt_lc))]
        time_windows2 = [[k - dt_lc, k] for k in
                         np.abs(np.arange(-args['mjd']-dt_lc/2, -MJD[0], dt_lc))]
        time_windows.extend(time_windows2) 

    if (time_windows[-1][1] - MJD[0]) < dt_lc:
        del time_windows[-1]
    time_windows.append('')

    for t_window in time_windows:
        if t_window == '':
            opath = src.sed_path
            opath_1GeV = opath+'_1GeV'
            partition='long'
        else:
            add_str = '{:.1f}_{:.1f}'.format(t_window[0], t_window[1])
            opath = os.path.join(src.lc_path, add_str)
            opath_1GeV = os.path.join(src.lc_path+'_1GeV', add_str)
            partition='kta'
        submit_fit(sub_args, opath, src, sub_file=ev_str+'.sub', trange=t_window, partition=partition, **args)
        if args['emin'] < 1e3:
            opath = os.path.join(src.lc_path, add_str)
            submit_fit(tsub_args, opath_1GeV, src, sub_file=ev_str+'.sub', trange=t_window, partition=partition, **args)
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
            plot.make_ts_plot(ts_map_path, os.path.join(vou_out, 'src_dict.npy'),
                              os.path.join(vou_out, 'find_out_temp.txt'),
                              plt_mode='tsmap', yaxis=yaxis, **args)
        except Exception as inst:
            warnings.warn("Couldn't create ts map...")
            print(inst)

        try:
            plot.make_ts_plot(ts_map_short_path, os.path.join(vou_out, 'src_dict.npy'),
                              os.path.join(vou_out, 'find_out_temp.txt'),
                              plt_mode='tsmap', legend=False, **args)
        except Exception as inst:
            warnings.warn("Could't create residual map...")
            print(inst)

        try:
            plot.make_ts_plot(ts_map_path, os.path.join(vou_out, 'src_dict.npy'),
                              os.path.join(vou_out, 'find_out_temp.txt'),
                              plt_mode='residmap', legend=False, **args)
        except Exception as inst:
            warnings.warn("Couldn't create residual map...")
            print(inst)

        for src in analysis.srcs:
            print('Make Plots for Source {}'.format(key))
            try:
                plot.make_lc_plot(src.lc_path, args['mjd'])
            except Exception as inst:
                warnings.warn("Couldn't create lightcurve for source {}".format(key))
                print(inst)
            try:
                plot.make_lc_plot(src.lc_path + '_1GeV', args['mjd'])
            except Exception as inst:
                warnings.warn("Couldn't create additional 1GeV lightcurve for source {}".format(key))
                print(inst)

            folders = [fold for fold in os.listdir(src.lc_path) if
                       os.path.isdir(os.path.join(src.lc_path,fold))]
            for fold in folders:
                try:
                    seds_list = [(os.path.join(src.lc_path, fold), 'k', 'red' , True, True, True),(src.sed_path, 'grey',
                                'grey', True, True, True)]
                    if os.path.exists(os.path.join(src.lc_path+'_1GeV', fold)):
                        seds_list.append((os.path.join(src.lc_path + '_1GeV', fold), 'k', 'blue' , True, True, False))
                    plot.make_sed_plot(seds_list, mw_data=src.mc_data_path, dec=src.dec)
                except Exception as inst:
                    warnings.warn("Couldn't create SED for source {}".format(key))
                    print(inst)
                    pass
            try:
                plot.make_sed_plot([(src.sed_path, 'grey', 'grey', True, True, True)],
                                   mw_data=src.mw_data_path, dec=src.dec)
            except Exception as inst:
                    warnings.warn("Couldn't create all year SED")
                    print(inst)
                    pass
            try:
                make_gif(src.lc_path)
            except Exception as inst:
                warnings.warn("Couldn't create an animated light curve {}".format(key))
                print(inst)
    src_latex = ''
    dists = np.array([s['dist'] for s in src_dict])
    sinds = np.argsort(dists)
    for src in src_dict[sinds]:
        if src['dist'] > args['max_dist']:
            print('Source exceeds distance of {} deg. No summary information will be added'.format(args['max_dist']))
            continue
        src_latex += source_summary(bpath, src, args['mjd'], mode=args['mode'])
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

    c = SkyCoord(args['ra'], args['dec'], frame='icrs', unit="deg")
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
                          ra=args['ra'],
                          dec=args['dec'],
                          emin=1.*args['emin']/1000.,
                          l=gal.l.deg,
                          b=gal.b.deg,
                          cat_srcs=short_out,
                          rx_map=os.path.join(vou_out, 'RX_map.eps'),
                          vou_pic=os.path.join(vou_out, 'candidates.eps'),
                          ts_map=os.path.join(bpath, 'ts_map/tsmap.png'),
                          ts_map_short=os.path.join(bpath, 'ts_map_short/tsmap.png'),
                          vou_output=full_out,
                          event=ev_str,
                          tsmjd1=tsmjd1,
                          tsmjd2=tsmjd2,
                          src_latex=src_latex,
                          mjd1=MJD[0],
                          mjd2=MJD[1],
                          energy=args['emin']/1000.,
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
