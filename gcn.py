import numpy as np
import requests
from astropy.time import Time


def parse_direction(string):
    if isinstance(string, float):
        return string
    else:
        return float(string.split('{')[0].replace('d',''))
    
def parse_error(string):
    if isinstance(string, float):
        return string
    else:
        return float(string.split('[')[0])/60.

def make_mjd(string1, string2):
    s1_split = string1.split(';')[2].split('(')[0].replace(' ', '').replace('/', '-')
    s2_split = string2.split('{')[1].split('}')[0].replace(' ', '')
    t = Time('20'+s1_split+' '+s2_split, format='iso', scale='utc')
    return t.mjd

def calc_signess(string):
    if isinstance(string, float):
        return string
    else:
        return float(string.split('[')[0])

def read_gcn(url, basepath=None):
    split_sign = '//////////////////////////////////////////////////////////////////////'
    html= requests.get(url).text.encode('utf8')
    html_split = html.split(split_sign)
    for c in range(1,len(html_split)):
        split0 = html_split[c].split('\n')
        split_info = []
        for i in split0:
            i_split = i.split(':   ')
            if len(i_split) > 1:
                split_info.append([i_split[0].replace(':',''), i_split[-1].strip()])
        keys = ['NOTICE_DATE', 'NOTICE_TYPE', 'SRC_RA', 'SRC_DEC', 'SRC_ERROR', 'SRC_ERROR50', 'DISCOVERY_DATE',
                'SIGNALNESS', 'SUN_DIST', 'MOON_DIST', 'DISCOVERY_TIME']
        gcn_info = dict()
        for i in split_info:
            if i[0] in keys:
                gcn_info[i[0]] = i[1] 
        gcn_info['SRC_DEC'] = parse_direction(gcn_info['SRC_DEC'])
        gcn_info['SRC_RA'] = parse_direction(gcn_info['SRC_RA'])
        gcn_info['SRC_ERROR'] = parse_error(gcn_info['SRC_ERROR'])
        gcn_info['SRC_ERROR50'] = parse_error(gcn_info['SRC_ERROR50'])
        gcn_info['SIGNALNESS'] = calc_signess(gcn_info['SIGNALNESS'])
        gcn_info['MJD'] = make_mjd(gcn_info['DISCOVERY_DATE'], gcn_info['DISCOVERY_TIME'])
        if basepath is not None:
            np.save(os.path.join(basepath, 'gcn_info_{}.npy'.format(c)), gcn_info)
    return gcn_info
