# coding: utf-8

import sys
sys.path.append('/scratch9/tglauch/Fermi_Tools/get_fermi_data')
sys.path.append('/scratch9/tglauch/realtime_service/python/lib/python2.7/site-packages')
sys.path.append('/home/ga53lag/Software/python_scripts/')
from get_fermi_data import get_data
import subprocess
import argparse
import os
from roi_functions import get_sources
import shutil
import datetime
from slack_lib import print_to_slack
import numpy as np
from myfunctions import MET_to_MJD, dict_to_nparray
import plot
import warnings
import re
import time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time

path_settings = {'sed': 'all_year/sed',
                 'lc': 'lightcurve'}

partition_t = {'kta':'2:29:59', 'long':'1-23:59:59'}



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
        type=float, default=28)
    parser.add_argument(
        "--dt",
        help="length of time window to analyze",
        type=float, default=28)
    parser.add_argument(
        "--emin",
        help="lower energy bound for SED",
        type=float, default=1000)
    parser.add_argument(
        "--event",
        help="event name",
        required=False)
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
    args = parser.parse_args()
    return args.__dict__


def source_summary(bpath, src):
    print src['name']
    bpath_src = src_path(bpath, src['name'])
    lc_path = os.path.join(bpath_src, path_settings['lc'], 'lightcurve.pdf')
    sed_path = os.path.join(bpath_src, path_settings['sed'])
    l_str ='\subsection{{{srcname}}}\n'.format(srcname=src['name'])
    try:
        fit_res = np.load(os.path.join(sed_path, 'llh.npy'))[()]
        sed_path = os.path.join(sed_path, 'sed.pdf')
        ts = fit_res['sources'][src['name']]['ts']
        ra = fit_res['sources'][src['name']]['RAJ2000']
        dec = fit_res['sources'][src['name']]['DEJ2000']
        t_str = 'ra = {ra:.2f}$^\circ$ , dec = {dec:.2f}$^\circ$ , ts = {ts:.2f}, distance = {dist:.2f}$^\circ$ \n'
        l_str += t_str.format(ts=ts, ra=ra, dec=dec, dist=src['dist'])
    except Exception as inst:
        warnings.warn("Can not find fit result for {}".format(src['name']))
        print(inst)

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
            prob_str += '{:.2f} & {:.2f} & {:.2f} \\\\ \n'.format(mjd[prob_mask][ind][i],
                                                                  src_prob[prob_mask][i]*100,
                                                                  energy[prob_mask][i]/1e3)
        tab_str = tab_str.format(events=prob_str)
        l_str += tab_str
    except Exception as inst:
        warnings.warn("Could not find source probabilities")
        print(inst)
    print l_str
    with open('source_summary.tex', 'r') as infile:
        fig_str = infile.read()       
    if os.path.exists(sed_path):
        l_str += fig_str.format(path = sed_path, caption='SED for {}'.format(src['name']))
    if os.path.exists(lc_path):
        l_str += fig_str.format(path = lc_path, caption='Lightcurve for {}'.format(src['name']))
    l_str += '\\clearpage \n'
    return l_str


def submit_fit(args, opath, src_arr, trange='', ana_type='SED', partition='kta', **kwargs):
    if kwargs.get('make_pdf'):
        return    
    print('submit with args {}'.format(args))
    odtype = np.dtype([('name', np.unicode, 32), ('ra', np.float32), ('dec', np.float32)])
    if isinstance(src_arr, dict):
        src_arr=dict_to_nparray(src_arr, dtype=odtype)
    src_arr = np.atleast_1d(src_arr)
    if trange != '':
        args += ' --time_range {} {} '.format(trange[0], trange[1])
    args += ' --outfolder {} '.format(opath)
    if not os.path.exists(opath):
        os.makedirs(opath)
    mask = np.array(['3FGL' not in sname for sname in src_arr['name']])
    if len(src_arr[mask])>0:
        xml_path = os.path.join(opath, 'add_source.xml')
        generate_src_xml(src_arr[mask], xml_path)
        args += '--xml_path {}'.format(xml_path)
    with open('./slurm_draft.xml', 'r') as f:
        submit = (f.read()).format(bpath=opath, args=args,
                                   ana_type=ana_type,
                                   time=partition_t[partition],
                                   partition=partition)
    submitfile = os.path.join(opath, 'fermi.sub')
    with open(submitfile, "w+") as file:
        file.write(submit)
    os.system("sbatch {}".format(submitfile))
    return opath


def generate_src_xml(src_arr, xml_path):
    xml_str =\
'<?xml version="1.0" ?>\n\
<source_library title="source library">\n\
<!-- Point Sources -->\n'
    with open('./src_xml_draft.xml', 'r') as f:
        xml_temp = f.read()
    for src in src_arr:
        print('Generate Source {}'.format(src['name']))
        xml_str += xml_temp.format(ra=src['ra'], dec=src['dec'], name=src['name'])
        xml_str += '\n'
    xml_str+='</source_library>'
    with open(xml_path, 'w+') as f:
        f.write(xml_str)
    return


def src_path(bpath, src):
    return os.path.join(bpath, src.replace(' ', '_'))


args = parseArguments()
rec = args['recovery']
make_pdf = args['make_pdf']
print('Run with args')
print(args)
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
bpath = '/scratch9/tglauch/realtime_service/output/{}'.format(ev_str)
this_path = os.path.dirname(os.path.abspath(__file__))
fermi_data = os.path.join(bpath, 'fermi_data')
vou_out = os.path.join(bpath, 'vou_blazar')

if not rec and not make_pdf:
    cmd = [os.path.realpath('/scratch9/tglauch/VOU_Blazars/bin/vou-blazars'),
           str(args['ra']), str(args['dec']), str(120), str(30), str(90)]
    if not os.path.exists(vou_out):
        os.makedirs(vou_out)
    os.chdir(vou_out)
    subprocess.call(cmd)

    # Setup Variables
    src_dict, out_str = get_sources(args['ra'], args['dec'])
    np.save('src_dict.npy', src_dict)
    print_to_slack(out_str)

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
        print_to_slack(lines, os.path.join(vou_out, 'candidates.eps'))
    with open('./full_output', 'w+') as ofile:
        ofile.write(lines.replace('\n' ,'\\\ \n'))
    os.chdir(this_path)

    # download gamma-ray data
    kwargs = {'emin': args['emin'], 'days': args['dt'], 'out_dir' : fermi_data}
    if args['mjd_range'] is None:
        kwargs['mjd'] = args['mjd']
        kwargs['mode'] = args['mode']
    else:
        kwargs['mjd'] = args['mjd_range']
    print kwargs
    MET = get_data(args['ra'], args['dec'], **kwargs)
    MJD = [MET_to_MJD(float(i)) for i in MET]
    run_info = {'MJD' : MJD,
                'args' : args,
                'src_dict' : src_dict}
    np.save(os.path.join(bpath,'run_info.npy'), run_info)
else:
    run_info = np.load(os.path.join(bpath,'run_info.npy'))[()]
    MJD = run_info['MJD']
    args = run_info['args']
    args['recovery'] = rec
    args['make_pdf'] = make_pdf
    src_dict = run_info['src_dict']
print src_dict

# start the gamma-ray analysis

#TS maps
sargs = '--target_src center --free_radius 2 --data_path {} --use_3FGL --emin {} '
sargs = sargs.format(fermi_data, args['emin'])
ts_dict = {'name': ['center'], 'ra': [args['ra']], 'dec': [args['dec']]}
ts_map_path = os.path.join(bpath, 'ts_map')
submit_fit(sargs, ts_map_path, ts_dict, ana_type='TS_Map', partition='long', **args)

#getsrcprob
sargs = '--target_src center --free_radius 2 --data_path {} --use_3FGL --emin {} '
sargs = sargs.format(fermi_data, args['emin'])
ts_dict = {'name': np.concatenate([['center'],src_dict['name']]), 
           'ra': np.concatenate([[args['ra']], src_dict['ra']]),
           'dec': np.concatenate([[args['dec']], src_dict['dec']])
          }
srcprob_path = os.path.join(bpath, 'srcprob')
submit_fit(sargs, srcprob_path, ts_dict, ana_type='srcprob', partition='long', **args)

print('Submit_SEDs')
job_dict = {}
for src in src_dict:
    bpath_src = src_path(bpath, src['name'])
    print bpath_src
    if os.path.exists(bpath_src) and (not args['recovery']) and (not args['make_pdf']):
        print('Remove Path: {}'.format(bpath_src))
        shutil.rmtree(bpath_src)
    sargs = '--target_src {} --free_radius 2 --data_path {} --use_3FGL --emin {} '
    sargs = sargs.format(src['name'].replace(' ', '_'), fermi_data, args['emin'])
    print args
    time_windows = [[k, k + args['dt_lc']] for k in
                    np.arange(MJD[0], MJD[1], args['dt_lc'])]
    time_windows.append('')
    job_dict[src['name']] = {'sed': os.path.join(bpath_src, path_settings['sed']),
                             'lc': os.path.join(bpath_src, path_settings['lc'])}
    print job_dict
    for t_window in time_windows:
        print('{}'.format(t_window))
        if t_window == '':
            opath = job_dict[src['name']]['sed']
            partition='long'
        else:
            add_str = '{:.1f}_{:.1f}'.format(t_window[0], t_window[1])
            opath = os.path.join(job_dict[src['name']]['lc'], add_str)
            partition='kta'
        submit_fit(sargs, opath, src, trange=t_window, partition=partition, **args)

mins = 0
final_pdf = False
while not final_pdf:
    if not make_pdf:
        time.sleep(60)
    mins += 1
    jobs = os.popen('squeue --user ga53lag').read()
    len_jobs = len([i for i in jobs.split('\n') if 'fermi' in i])
    print len_jobs
    if len_jobs == 0 or make_pdf:
        final_pdf = True
    if not mins % 15 == 1 and not final_pdf:
        continue

    try:
        plot.make_ts_plot(ts_map_path, os.path.join(vou_out, 'src_dict.npy'),
                          os.path.join(vou_out, 'candidates_image_position.txt'),
                          mode='tsmap')
    except Exception as inst:
        warnings.warn("Couldn't create ts map...")
        print(inst)

    try:
        plot.make_ts_plot(ts_map_path, os.path.join(vou_out, 'src_dict.npy'),
                          os.path.join(vou_out, 'candidates_image_position.txt'),
                          mode='residmap')
    except Exception as inst:
        warnings.warn("Couldn't create residual map...")
        print(inst)

    for key in job_dict.keys():
        print('Make Plots for Source {}'.format(key))
        try:
            plot.make_lc_plot(job_dict[key]['lc'])
        except Exception as inst:
            warnings.warn("Couldn't create lightcurve for source {}".format(key))
            print(inst)
        try:
            plot.make_sed_plot(job_dict[key]['sed'])
        except Exception as inst:
            warnings.warn("Couldn't create SED for source {}".format(key))
            print(inst)

    src_latex = ''
    for src in src_dict:
        src_latex += source_summary(bpath, src)
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
                          res_map=os.path.join(bpath, 'ts_map/residmap.png'),
                          vou_output=full_out,
                          event=ev_str,
                          src_latex=src_latex,
                          mjd1=MJD[0],
                          mjd2=MJD[1],
                          energy=args['emin']/1000.)
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
