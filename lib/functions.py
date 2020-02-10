# coding: utf-8

import sys
sys.path.append('/home/ga53lag/Software/python_scripts/')
import numpy as np
from astropy.coordinates import SkyCoord
import os
import collections
from myfunctions import dict_to_nparray
from astropy.io import fits
import imageio
from scipy.interpolate import interp1d

def get_symbol(src_code):
    src_code = str(src_code)
    if src_code == '-1111':
        return False # gamma sources should be cataloged in Fermi so no plotting here
        return dict(marker = "^", color='blue', zorder=25,
                    ms=7, fillstyle='none')
    elif (src_code == '-99999') or (src_code == '-88888'):
        return dict(marker = "p", color='purple')
    else:
        if src_code[0] =='-':
            return False
            mcode = src_code[0:2]
        else:
            return False
        return cat_dict[int(mcode)][1]
    
def radio_circle(src_code):
    src_code = str(src_code)
    if (src_code == '-1111') or (src_code == '-99999') or (src_code == '-88888'):
        return False
    if src_code[0] =='-':
        mcode = src_code[0:2]
        if np.abs(int(mcode)) not in [3,7]:
            return False
    else:
        mcode = src_code[0]
    if cat_dict[int(mcode)][2] is not False:
        return cat_dict[int(mcode)][2]
    else:
        return False
    

def xray_circle(src_code):
    src_code = str(src_code)
    if (src_code == '-1111') or (src_code == '-99999') or (src_code == '-88888'):
        return False
    if src_code[0] =='-':
        return False
        mcode = src_code[0:2]
    else:
        mcode = src_code[0]
    if cat_dict[int(mcode)][3] is not False:
        return cat_dict[int(mcode)][3]
    else:
        return False
    return
    
def get_circle_size(code):
    om = int(code/10000.)
    code = abs(code)
    cxray  = int((code-abs(om)*10000.)/100.)
    cradio = code-abs(om)*10000.-cxray*100.
    if cradio > 0:
        cs = max(1.0,cradio*8./99.)
    else:
        cs = 0
    if cxray > 0:
        cx = max(1.0,cxray*8./99.)
    else:
        cx =0 
    return cs, cx



cat_dict = { 1: ['HSP', dict(marker = ".", zorder=10), 'orange', 'orange'],
            2: ['ISP', dict(marker = ".", zorder=7), 'cyan', 'darkcyan'],
            3: ['LSP',dict(marker = ".", zorder=3), 'blue', 'purple'],
            4: ['non-jet AGN',dict(marker = ".", color='red', ms =0), False, False] ,
            5:['unkown type', dict(marker = ".", color='red', ms =0), False, False],
            -1 :['3HSP', dict(marker = "*", color='grey',  ms=10, markeredgecolor='k', zorder=15),
                 False, False],
            -2 : ['BZCat', dict(marker = "d", color='grey', ms=10, markeredgecolor='k', zorder=20),
                  False, False],
            -3 : ['Crates', dict(marker = "s", color='grey', fillstyle='none', zorder=5, ms=7),
                  'grey', False],
            -4: ['Cluster', dict(marker = ".", color='red', ms =3, zorder=5),
                 False, False],
            -5 : ['3HSP w/o match',
                  dict(marker = "*", fillstyle='none', color='grey', ms=10, zorder=15),
                  False, False],
            -6: ['5BZCat w/o match', dict(marker = "d", ms=10, color='grey', zorder=20),
                 False, False],
            -7: ['Crates w/o match', dict(marker = "s", color='blue', fillstyle='none', zorder=5, ms=7),
                 'grey', False],
            -8: ['X-ray sources', dict(marker = ".", color='red', ms =0), False, False,],
            -9: ['radio sources', dict(marker = ".", color='red', ms =0), False, False],
            -10: ['fermi sources', dict(marker = "^", color='grey', zorder=25, ms=10, fillstyle='none'), False, False]}


 
fields = ['name', 'ra', 'dec', 'alt_name']
odtype = np.dtype([('name', np.unicode, 32), ('ra', np.float32), ('dec', np.float32), ('dist', np.float32)]) 
path_settings = {'sed': 'all_year/sed',
                 'lc': 'lightcurve'}
#vou_path = '/scratch9/tglauch/Software/VOU_Blazars/v2/bin/vou-blazars'
vou_path ='/home/ga53lag/Software/VOU_Blazars/v5/VOU_Blazars//bin/vou-blazars'
partition_t = {'kta':'2:30:00', 'long':'2-00:00:00', 'xtralong': '7-00:00:00'}

def get_68_psf(E):
    x = np.genfromtxt('./lib/lat_68_psf.txt', delimiter = ',')
    return interp1d(x[:,0], x[:,1])(E)


def submit_fit(args, opath, srcs=None, trange='', sub_file='fermi.sub',  ana_type='SED', partition='kta'):
    if trange != '':
        args += ' --time_range {} {} '.format(trange[0], trange[1])
    if not os.path.exists(opath):
        os.makedirs(opath)
    args += ' --outfolder {} '.format(opath)
    if srcs is not None:
        non_fermi_srcs = [src for src in srcs if not '4FGL' in src.name]
        if len(non_fermi_srcs)>0:
            xml_path = os.path.join(opath, 'add_source.xml')
            generate_src_xml(non_fermi_srcs, xml_path)
            args += '--xml_path {}'.format(xml_path)
    with open('./lib/slurm_draft.xml', 'r') as f:
        submit = (f.read()).format(bpath=opath, args=args,
                                   ana_type=ana_type,
                                   time=partition_t[partition],
                                   partition=partition)
    submitfile = os.path.join(opath, sub_file)
    with open(submitfile, "w+") as file:
        file.write(submit)
    print('submit with args {}'.format(args))
    os.system("sbatch {}".format(submitfile))
    return opath


def generate_src_xml(srcs, xml_path):
    xml_str =\
'<?xml version="1.0" ?>\n\
<source_library title="source library">\n\
<!-- Point Sources -->\n'
    with open('./lib/src_xml_draft.xml', 'r') as f:
        xml_temp = f.read()
    for src in srcs:
        xml_str += xml_temp.format(ra=src.ra, dec=src.dec, name=src.name)
        xml_str += '\n'
    xml_str+='</source_library>'
    with open(xml_path, 'w+') as f:
        f.write(xml_str)
    return


def get_lc_time3fgl(src_of_interest, emin=1e3):
    data = fits.open('./lib/gll_psc_v16.fit')
    src_ind = np.where(data[1].data['Source_Name']==src_of_interest)[0][0]
    
    E_bins = [[100, '100_300'],
              [300, '300_1000'],
              [1000, '1000_3000'], 
              [3000, '3000_10000'], 
              [10000, '10000_100000']]
    
    for i in range(len(E_bins)):
        if E_bins[i][0] >= emin:
            E_bins[i][0] = 1
        elif E_bins[i][0] < emin and E_bins[i+1][0] > emin:
            E_bins[i][0] = 1.*(E_bins[i+1][0]-emin)/(E_bins[i+1][0] - E_bins[i][0])
        else:
            E_bins[i][0] = 0
    ts_sum = np.sum([(data[1].data[src_ind]['Sqrt_TS{}'.format(e_rang[1])])**2*e_rang[0] 
                     for e_rang in E_bins])
    t_sens = (6/ts_sum)*(365*4-4)
    print t_sens
    if t_sens < 28.:
        t_disc = (36/ts_sum)*(365*4-4)
        print t_disc
        if t_disc < 28.:
            return t_disc
        else:
            return 28.
    return t_sens

# Thanks to the improved statistics, the source photon fluxes in 4FGL are reported in seven energybands (1: 50 to 100
# MeV; 2: 100 to 300 MeV; 3: 300 MeV to 1 GeV; 4: 1 to 3 GeV; 5: 3 to 10 GeV; 6:10 to 30 GeV; 7:  30 to 300 GeV)

def get_lc_time4fgl(src_of_interest, emin=1e3):
    data = fits.open('./lib/gll_psc_v19.fit')
    src_ind = np.where(data[1].data['Source_Name']==src_of_interest)[0][0]

    E_bins = [50,100,300,1000,3000,10000,30000,300000]

    for i in range(len(E_bins)):
        if E_bins[i] >= emin:
            E_bins[i] = 1
        elif E_bins[i] < emin and E_bins[i+1] > emin:
            E_bins[i] = 1.*(E_bins[i+1]-emin)/(E_bins[i+1] - E_bins[i])
        else:
            E_bins[i] = 0
    ts_sum = np.sum([data[1].data[src_ind]['Sqrt_TS_Band'][i]**2*E_bins[i]
                     for i in range(len(E_bins)-1)])
    t_sens = (6/ts_sum)*(365*8)
    if t_sens < 100:
        t_disc = (25/ts_sum)*(365*8)
        t_disc = np.max([t_disc, 56])
        t_disc = np.min([t_disc, 100])
        return t_disc
    print(t_sens)
    return np.min([t_sens, 365])

def GreatCircleDistance(ra_1, dec_1, ra_2, dec_2):
    '''Compute the great circle distance between two events'''
    '''All coordinates must be given in radians'''
    delta_dec = np.abs(dec_1 - dec_2)
    delta_ra = np.abs(ra_1 - ra_2)
    x = (np.sin(delta_dec / 2.))**2. + np.cos(dec_1) *\
        np.cos(dec_2) * (np.sin(delta_ra / 2.))**2.
    return 2. * np.arcsin(np.sqrt(x))


def make_gif(basepath):
    images = []
    filenames = sorted([os.path.join(basepath, i, 'sed.pdf') for i in os.listdir(basepath) if
                        os.path.isdir(os.path.join(basepath, i))])
    opath = os.path.join(basepath,'movie.gif')
    if os.path.exists(opath):
        os.remove(opath)
    for filename in filenames:
        if not os.path.exists(filename):
            continue
        filename_png = filename.replace('.pdf', '.png')
        if os.path.exists(filename_png):
            os.remove(filename_png)
        os.system('convert -density 150 {} -quality 90 {}'.format(filename, filename_png))
        images.append(imageio.imread(filename.replace('.pdf', '.png')))
    if len(images)>0:
        imageio.mimsave(opath, images, duration=0.75)
    return


