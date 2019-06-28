import numpy as np
import os
import argparse
from roi_functions import submit_fit, generate_src_xml, path_settings, vou_path, partition_t
from myfunctions import GreatCircleDistance

def parseArguments():
    """Parse the command line arguments
    Returns:
    args : Dictionary containing the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ra",
        help="right ascension",
        type=float,
        required=True)
    parser.add_argument(
        "--dec",
        help="declination",
        type=float,
        required=True)
    parser.add_argument(
        "--event",
        help="event name",
        required=True)
    parser.add_argument(
        "--vou_name",
        help="preliminary vou name",
        required=True)
    args = parser.parse_args()
    return args.__dict__


nargs = parseArguments()
bpath= '/scratch9/tglauch/realtime_service/output/{}'.format(nargs['event'])
fermi_data = os.path.join(bpath, 'fermi_data')
run_info = np.load(os.path.join(bpath,'run_info.npy'), allow_pickle=True)[()]
args = run_info['args']
src = nargs['vou_name']
ev_str = args['event']
MJD = run_info['MJD']
src_folder = os.path.join(bpath, src)
sed_path = os.path.join(bpath, src,  path_settings['sed'])
lc_path = os.path.join(bpath, src, path_settings['lc'])
if not os.path.exists(src_folder):
    os.makedirs(src_folder)
if not os.path.exists(lc_path):
    os.makedirs(lc_path)
if not os.path.exists(sed_path):
    os.makedirs(sed_path)
dist = GreatCircleDistance(nargs['ra'], nargs['dec'], args['ra'],
                           args['dec'], unit='deg')
dtype=run_info['src_dict'].dtype
src_arr = np.array([(src, nargs['ra'], nargs['dec'], dist)], dtype=dtype)
run_info['src_dict'] = np.append(run_info['src_dict'][:2],  src_arr)
print('Source is at a distance of {}'.format(dist))
np.save(os.path.join(bpath,'run_info.npy'), run_info)
os.system('{vou_path} {ra} {dec} {loc_str} -s ; cat Sed.txt > {bpath}'.format(vou_path=vou_path, ra=nargs['ra'], dec=nargs['dec'],
                                                                              bpath=os.path.join(src_folder, 'sed.txt'), loc_str=2))
sargs = '--target_src {} --free_radius 2 --data_path {} --use_4FGL --emin {} '
sub_args = sargs.format(src.replace(' ', '_'), fermi_data, args['emin'])
tsub_args = sargs.format(src.replace(' ', '_'), fermi_data, 1e3)
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
        opath = sed_path
        opath_1GeV = opath+'_1GeV'
        partition='long'
    else:
        add_str = '{:.1f}_{:.1f}'.format(t_window[0], t_window[1])
        opath = os.path.join(lc_path, add_str)
        opath_1GeV = os.path.join(lc_path+'_1GeV', add_str)
        partition='kta'
    submit_fit(sub_args, opath, src_arr, sub_file=ev_str+'.sub', trange=t_window, partition=partition, **args)
    if args['emin'] < 1e3:
        opath = os.path.join(lc_path, add_str)
        submit_fit(tsub_args, opath_1GeV, src_arr, sub_file=ev_str+'.sub', trange=t_window, partition=partition, **args)


