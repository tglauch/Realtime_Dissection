from roi_functions import GreatCircleDistance, get_add_info, path_settings
import numpy as np
from source_class import Source 
import collections
import os
from astropy.coordinates import SkyCoord

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

class Analysis(object):
    def __init__(self, ev_name, ra, dec):
        self.event_name = ev_name
        self.ra = ra  
        self.dec = dec
        self.srcs = []
        
    def make_pdf(self):
        # To be implemented
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
    
