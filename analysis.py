from roi_functions import GreatCircleDistance, get_add_info, path_settings, vou_path
import numpy as np
from source_class import Source 
import collections
import os
from astropy.coordinates import SkyCoord
from gcn import read_gcn
from astropy.time import Time
from myfunctions import MET_to_MJD
from get_fermi_data import get_data
import subprocess
from slack_lib import print_to_slack
import re
import string
import random

files = collections.OrderedDict([
     ('4fgl', {'file': '4fgl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('3fgl', {'file': '3fgl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('3fhl', {'file': '3fhl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('5bzcat', {'file': '5bzcat.1.csv',
                'keys': ['Name', 'RAJ2000', 'DEJ2000']}),
     ('3hsp', {'file': '3hsp.1.csv',
              'keys': ['Name', 'ra', 'dec']}),
     ('fermi8yr', {'file': 'fermi8yr.1.csv',
                  'keys': ['Source_Name', 'RAJ2000', 'DEJ2000']}),
     ('crates', {'file': 'crates.1.csv',
                'keys': ['name', 'ra', 'dec']})])

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class Analysis(object):
    def __init__(self, mjd=None, ra=None, dec=None):
        self.event_name = None
        self.ra = ra  
        self.dec = dec
        self.srcs = []
        self.err90 = 2.0
        self.mode = 'end'
        self.mjd = mjd
        self.mjd_range = None
        self.gcn = 'Not Given'
        self.vou_out = None
        self.fermi_data = None
        self.ts_maps = []
        self.srcprob_path = None
        self.id = id_generator(size=5) 

    def make_pdf(self):
        # To be implemented
        return

    def ROI_analysis(self, this_path, radius):
        # To be implemented
        self.vou_out = os.path.join(self.bpath, 'vou_blazar')
        if not os.path.exists(self.vou_out):
            os.makedirs(self.vou_out)
        os.chdir(self.vou_out)
        cmd = [vou_path,
               str(self.ra), str(self.dec), radius, str(30), str(90)]
        for i in range(2):
            subprocess.call(cmd)

        # Setup Variables
        out_str = self.get_sources()
        headline = '*Result for {} *\n'.format(self.event_name)
        sum_text = headline + out_str
        print_to_slack(sum_text)

        # Convert VOU Blazar Output
        rx_ps = os.path.join(self.vou_out, 'RX_map.ps')
        os.system('ps2eps -B ' + rx_ps)
        cand_ps = os.path.join(self.vou_out, 'candidates.ps')
        os.system('ps2eps -B ' + cand_ps )

        #Create VOU Source Summary
        with open('./short_output', 'w+') as ofile:
            ofile.write(out_str.replace('\n' ,'\\\ \n'))
        with open('./phase1', 'r') as ifile:
            lines = ifile.read().split('Gamma-ray Counterparts')[0]
            lines = re.sub('\\[..?;.?m', ' ', lines)
            lines = lines.replace('[0m', ' ')
            print_to_slack('', pic=os.path.join(self.vou_out, 'candidates.eps'))
        with open('./full_output', 'w+') as ofile:
            ofile.write(lines.replace('\n' ,'\\\ \n'))
        os.chdir(this_path)
        return

    def get_fermi_data(self, days=None, mjd_range=None, overwrite=False):
        if (self.fermi_data != None) & (overwrite is False):
            return ## Data is already there
        self.fermi_data = os.path.join(self.bpath, 'fermi_data')
        args = {'emin': self.emin,
                'out_dir' : self.fermi_data}
        if days is not None:
            args['days'] = days
        if mjd_range is None:
            args['mjd'] = self.mjd
            args['mode'] = self.mode
        else:
            args['mjd'] = mjd_range
        MET = get_data(self.ra, self.dec, **args)
        MJD = [MET_to_MJD(float(i)) for i in MET]
        self.mjd_range = MJD
        return

    def update_gcn(self):
        self.from_gcn(self.gcn)
        return

    def from_gcn(self, url):
        gcn_dict = read_gcn(url)
        self.ra = gcn_dict['SRC_RA']
        self.dec = gcn_dict['SRC_DEC']
        self.err90 = gcn_dict['SRC_ERROR']
        self.err50 = gcn_dict['SRC_ERROR50']
        self.mjd = gcn_dict['MJD']
        self.gcn = url
        t = Time.now()
        if np.abs(self.mjd - t.mjd)<2:
            self.mode = 'end'
        else:
            self.mode = 'mid'
        return 
        
    def calc_src_distances_to_bf(self):
        for src in self.srcs:
            src.dist = src.angular_distance(self.ra, self.dec)
        return
            
    def sort_sources_by_distance(self):
        dists = np.array([i.dist for i in self.srcs])
        inds = np.argsort(dists)
        self.srcs = [self.srcs[j] for j in inds]
        return
    
    def get_sources(self):
        for key in list(files):
            if not os.path.exists(files[key]['file']):
                continue
            temp = np.genfromtxt(files[key]['file'], dtype=None,
                                 names=True, delimiter=',')
            temp = np.atleast_1d(temp)
            if key == '5bzcat':
                for src in temp:
                    ra_split = src[1].split(' ')
                    dec_split = src[2].split(' ')
                    c = SkyCoord('{}h{}m{}s'.format(ra_split[0],
                                                    ra_split[1],
                                                    ra_split[2]),
                                 '{}d{}m{}s'.format(dec_split[0],
                                                    dec_split[1],
                                                    dec_split[2]),
                                 frame='icrs')

                    src[1] = c.ra.degree
                    src[2] = c.dec.degree
                dtype = [('Name', 'S15'), ('RAJ2000', np.float32),
                         ('DEJ2000', '<f8')]
                temp = temp.astype(dtype)
            for src in temp:
                src_obj = Source(src[files[key]['keys'][0]],
                                 src[files[key]['keys'][1]],
                                 src[files[key]['keys'][2]])
                self.srcs.append(src_obj)

        i = 0
        while i < len(self.srcs):
            ang_diff = [self.srcs[i].angular_distance_to_source(self.srcs[j])
                        for j in range(i + 1, len(self.srcs))]
            to_delete = []
            for j, diff in enumerate(ang_diff):
                if diff < 0.2:
                    if not self.srcs[i + j + 1].name in self.srcs[i].names \
                        and self.srcs[i + j + 1].name != self.srcs[i].name:
                            self.srcs[i].names.append(self.srcs[i + j + 1].name)
                    to_delete.append(i + j + 1)
            for k, j in enumerate(to_delete):
                del(self.srcs[j - k])
            i += 1


        fmt_str = '*{}* {}\n ra: {:.2f} deg |  dec: {:.2f} deg | distance: {:.2f} deg [ra: {:.2f} deg , dec: {:.2f} deg]'
        self.calc_src_distances_to_bf()
        self.sort_sources_by_distance()
        out_str = ''
        for src in self.srcs:
            add_info = get_add_info(src.name)
            ostr= fmt_str.format(src.name, add_info, src.ra, src.dec,
                                 src.dist, src.ra - self.ra, src.dec - self.dec)
            if len(src.names) > 0:
                ostr +=  '\n Associations: {}'.format(', '.join(src.names))
            print(ostr)
            out_str += ostr + '\n\n'
        return out_str
    
