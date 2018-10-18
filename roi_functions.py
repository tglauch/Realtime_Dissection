import numpy as np
from astropy.coordinates import SkyCoord
import os
import collections

fields = ['name', 'ra', 'dec', 'alt_name']

files = collections.OrderedDict(
    {'3fgl': {'file': '3fgl.1.csv',
              'keys': ['name', 'ra', 'dec']},
     '3fhl': {'file': '3fhl.1.csv',
              'keys': ['name', 'ra', 'dec']},
     '5bzcat': {'file': '5bzcat.1.csv',
                'keys': ['Name', 'RAJ2000', 'DEJ2000']},
     '3hsp': {'file': '3hsp.1.csv',
              'keys': ['Name', 'ra', 'dec']},
     'fermi8yr': {'file': 'fermi8yr.1.csv',
                  'keys': ['Source_Name', 'RAJ2000', 'DEJ2000']},
     'crates': {'file': 'crates.1.csv',
                'keys': ['name', 'ra', 'dec']}})


def GreatCircleDistance(ra_1, dec_1, ra_2, dec_2):
    '''Compute the great circle distance between two events'''
    '''All coordinates must be given in radians'''
    delta_dec = np.abs(dec_1 - dec_2)
    delta_ra = np.abs(ra_1 - ra_2)
    x = (np.sin(delta_dec / 2.))**2. + np.cos(dec_1) *\
        np.cos(dec_2) * (np.sin(delta_ra / 2.))**2.
    return 2. * np.arcsin(np.sqrt(x))


def get_sources(ra, dec):
    src_dict = {}
    for key in fields:
        src_dict[key] = []
    for key in files.keys():
        if not os.path.exists(files[key]['file']):
            continue
        temp = np.genfromtxt(files[key]['file'], dtype=None,
                             names=True, delimiter=',')
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
        for i in range(len(files[key]['keys'])):
            app =np.atleast_1d(temp[files[key]['keys'][i]])
            src_dict[fields[i]].extend(app)

    i = 0
    while i < len(src_dict['dec']):
        src_dict['alt_name'].append([])
        ang_diff = [np.degrees(
                    GreatCircleDistance(np.radians(src_dict['ra'][i]),
                                        np.radians(src_dict['dec'][i]),
                                        np.radians(src_dict['ra'][j]),
                                        np.radians(src_dict['dec'][j])))
                    for j in range(i + 1, len(src_dict['dec']))]
        to_delete = []
        for j, diff in enumerate(ang_diff):
            if diff < 0.1:
                if not src_dict['name'][i + j + 1] in src_dict['alt_name'][i] \
                    and src_dict['name'][i + j + 1] != src_dict['name'][i]:
                        src_dict['alt_name'][i].append(src_dict['name'][i + j + 1])
                to_delete.append(i + j + 1)
        for k, j in enumerate(to_delete):
            del(src_dict['name'][j - k])
            del(src_dict['ra'][j - k])
            del(src_dict['dec'][j - k])
        i += 1

    fmt_str = '*{}* \n ra: {:.2f} deg |  dec: {:.2f} deg | distance: {:.2f} deg [ra:{:.2f} , dec:{:.2f}]'
    gcd = [GreatCircleDistance(np.radians(src_dict['ra'][i]), np.radians(src_dict['dec'][i]), np.radians(ra),
           np.radians(dec)) for i in range(len(src_dict['name']))]
    src_dict['dist'] = np.degrees(gcd)
    inds = np.argsort(gcd)
    out_str = '' 
    for i in inds:
        ostr= fmt_str.format(src_dict['name'][i],
                             src_dict['ra'][i],
                             src_dict['dec'][i],
                             np.degrees(gcd[i]),
                             src_dict['ra'][i] - ra,
                             src_dict['dec'][i] - dec)
        if len(src_dict['alt_name'][i]) > 0:
            ostr +=  '\n Alt Names: {}'.format(', '.join(src_dict['alt_name'][i]))
        print(ostr)
        out_str += ostr + '\n\n'
    return src_dict, out_str
