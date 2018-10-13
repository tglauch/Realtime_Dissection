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
from myfunctions import MET_to_MJD
import time
import plot
import warnings
import re


path_settings = {'sed': 'all_year/sed',
                 'lc': 'lightcurve'}


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
        type=float, default=2000)
    args = parser.parse_args()
    return args.__dict__


def source_summary(bpath, src):
    l_str = '\subsection{{{srcname}}} \
             ra = ${ra:.2f}^\circ$ , dec = ${dec:.2f}^\circ$ , ts = {ts:.2f}'
    with open('source_summary.tex', 'r') as infile:
        fig_str = infile.read()
    bpath_src = src_path(bpath, src)
    lc_path = os.path.join(bpath_src, path_settings['lc'], 'lightcurve.pdf')
    sed_path = os.path.join(bpath_src, path_settings['sed'])
    fit_res = np.load(os.path.join(sed_path, 'llh.npy'))[()]
    sed_path = os.path.join(sed_path, 'sed.pdf')
    ts = fit_res['sources'][src]['ts']
    print ts
    ra = fit_res['sources'][src]['RAJ2000']
    dec = fit_res['sources'][src]['DEJ2000']
    ## include right ascension and declination
    l_str = l_str.format(srcname=src, ts=ts, ra=ra, dec=dec)
    if os.path.exists(sed_path):
        l_str += fig_str.format(path = sed_path, caption='SED for {}'.format(src))
    if os.path.exists(lc_path):
        l_str += fig_str.format(path = lc_path, caption='Lightcurve for {}'.format(src))
    l_str += '\\clearpage \n'
    return l_str


def submit_fit(args, opath, src_dict, which, trange='', ana_type='SED'):
    if trange != '':
        args += ' --time_range {} {} '.format(trange[0], trange[1])
    args += ' --outfolder {} '.format(opath)
    if not os.path.exists(opath):
        os.makedirs(opath)
    if '3FGL' not in src_dict['name'][which]:
        xml_path = os.path.join(opath, 'add_source.xml')
        generate_src_xml(src_dict['name'][which], src_dict['ra'][which],
                         src_dict['dec'][which], xml_path)
        args += '--xml_path {}'.format(xml_path)
    with open('./slurm_draft.xml', 'r') as f:
        submit = (f.read()).format(bpath=opath, args=args, ana_type=ana_type)
    submitfile = os.path.join(opath, 'fermi.sub')
    with open(submitfile, "w+") as file:
        file.write(submit)
    os.system("sbatch {}".format(submitfile))
    return opath


def generate_src_xml(name, ra, dec, xml_path):
    with open('./src_xml_draft.xml', 'r') as f:
        xml_str = (f.read()).format(name=name, ra=ra, dec=dec)
    with open(xml_path, 'w+') as f:
        f.write(xml_str)
    return


def src_path(bpath, src):
    return os.path.join(bpath, src.replace(' ', '_'))


dtime = datetime.datetime.now()
ev_str = 'IC{}{:02d}{:02d}'.format(str(dtime.year)[-2:],
                                   dtime.month, dtime.day)
bpath = '/scratch9/tglauch/realtime_service/output/{}'.format(ev_str)
this_path = os.path.dirname(os.path.abspath(__file__))
fermi_data = os.path.join(bpath, 'fermi_data')

args = parseArguments()
cmd = [os.path.realpath('/scratch9/tglauch/VOU_Blazars/bin/vou-blazars'),
       str(args['ra']), str(args['dec']), str(90), str(30), str(60)]
vou_out = os.path.join(bpath, 'vou_blazar')
if not os.path.exists(vou_out):
    os.makedirs(vou_out)
os.chdir(vou_out)
subprocess.call(cmd)

print('Get Sources....')
src_dict, out_str = get_sources(args['ra'], args['dec'])
np.save('src_dict.npy', src_dict)
cand_ps = os.path.join(vou_out, 'candidates.ps')
cand_png = os.path.join(vou_out, 'candidates.png')
os.system('convert ' + cand_ps + ' -density 600 ' + cand_png)
print_to_slack(out_str)
with open('./phase1', 'r') as ifile:
    lines = ifile.read().split('Gamma-ray Counterparts')[0]
    lines = re.sub('\\[..?;.?m', ' ', lines)
    lines = lines.replace('[0m', ' ')
    print_to_slack(lines, cand_png)
os.chdir(this_path)
MET = get_data(args['ra'], args['dec'], emin=args['emin'],
               dt=args['dt'], out_dir=fermi_data)  # dt hardcoded!!!!
MJD = [MET_to_MJD(float(i)) for i in MET]
print MJD
print('Submit TS Map')
sargs = '--target_src center --free_radius 2 --data_path {} --use_3FGL --emin {} '
sargs = sargs.format(fermi_data, args['emin'])
ts_dict = {'name': ['center'], 'ra': [args['ra']], 'dec': [args['dec']]}
ts_map_path = os.path.join(bpath, 'ts_map')
submit_fit(sargs, ts_map_path, ts_dict, 0, ana_type='TS_Map')

print('Submit_SEDs')
job_dict = {}
for i, src in enumerate(src_dict['name']):
    bpath_src = src_path(bpath, src)
    print bpath_src
    if os.path.exists(bpath_src):
        shutil.rmtree(bpath_src)
    sargs = '--target_src {} --free_radius 2 --data_path {} --use_3FGL --emin {} '
    sargs = sargs.format(src.replace(' ', '_'), fermi_data, args['emin'])
    print args
    time_windows = [[k, k + args['dt_lc']] for k in
                    np.arange(MJD[0], MJD[1], args['dt_lc'])]
    time_windows.append('')
    job_dict[src] = {'sed': os.path.join(bpath_src, path_settings['sed']),
                     'lc': os.path.join(bpath_src, path_settings['lc'])}
    print job_dict
    for t_window in time_windows:
        print('{}'.format(t_window))
        if t_window == '':
            opath = job_dict[src]['sed']
        else:
            add_str = '{:.1f}_{:.1f}'.format(t_window[0], t_window[1])
            opath = os.path.join(job_dict[src]['lc'], add_str)
        submit_fit(sargs, opath, src_dict, i, trange=t_window)

len_jobs = 5
mins = 0
while (len_jobs > 0) and (mins < 45):
    print len_jobs
    jobs = os.popen('squeue --user ga53lag').read()
    len_jobs = len([i for i in jobs.split('\n') if 'fermi' in i])
    mins += 1
    time.sleep(60)

try:
    plot.make_ts_plot(ts_map_path, os.path.join(vou_out, 'src_dict.npy'),
                      mode='ts')
except Exception as inst:
    warnings.warn("Couldn't create ts map...")
    print(inst)

try:
    plot.make_ts_plot(ts_map_path, os.path.join(vou_out, 'src_dict.npy'),
                      mode='resmap')
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
for src in src_dict['name']:
    src_latex += source_summary(bpath, src)
with open(os.path.join(bpath, 'vou_blazar/full_output'), 'r') as f:
    full_out = f.read().encode('utf8').\
        replace('\x1b', '').replace('nu_p', 'nu peak')
    full_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', full_out)
    full_out = re.sub('(Match[^\\\\]*)', r'\\textbf{\1}', full_out)
    full_out = re.sub('(Dist[^\\\\]*)', r'\\textbf{\1}', full_out)
    full_out = [i.strip() + ' \n' for i in full_out.split('\n')]
    full_out = [i if i != '\\\\ \n' else '\\\\\\\ \n' for i in full_out]
    full_out = ' '.join(full_ out)
with open(os.path.join(bpath, 'vou_blazar/short_output'), 'r') as f:
    short_out = f.read().encode('utf8').replace('\x1b', '')
    short_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', short_out)
with open('./template.tex') as f:
    template = f.read()

out = template.format(cat_srcs=short_out,
                      vou_pic=os.path.join(bpath, 'vou_blazar/candidates.png'),
                      ts_map=os.path.join(bpath, 'ts_map/ts_map.png'),
                      res_map=os.path.join(bpath, 'ts_map/resmap_map.png'),
                      vou_output=full_out,
                      event=ev_str,
                      src_latex=src_latex)
with open('./summary.tex', 'w+') as f:
    f.write(out)

cmd = ['pdflatex', '-interaction', 'nonstopmode', './summary.tex']
proc = subprocess.Popen(cmd)
proc.communicate()

retcode = proc.returncode
if not retcode == 0:
    raise ValueError('Error {} executing command: {}'.format(retcode, ' '.join(cmd))) 