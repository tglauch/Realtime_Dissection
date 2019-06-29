# coding: utf-8

import sys
sys.path.append('/home/ga53lag/Software/python_scripts/')
import numpy as np
from astropy.coordinates import SkyCoord
import os
import collections
from myfunctions import dict_to_nparray
import scipy.integrate
import pyfits as fits 
import imageio

 
fields = ['name', 'ra', 'dec', 'alt_name']
odtype = np.dtype([('name', np.unicode, 32), ('ra', np.float32), ('dec', np.float32), ('dist', np.float32)]) 
path_settings = {'sed': 'all_year/sed',
                 'lc': 'lightcurve'}
vou_path = '/scratch9/tglauch/Software/VOU_Blazars/v2/bin/vou-blazars'
partition_t = {'kta':'2:30:00', 'long':'2-00:00:00', 'xtralong': '7-00:00:00'}

def submit_fit(args, opath, src_arr=None, trange='', sub_file='fermi.sub',  ana_type='SED', partition='kta', **kwargs):
    if kwargs.get('make_pdf'):
        return
    if trange != '':
        args += ' --time_range {} {} '.format(trange[0], trange[1])
    if not os.path.exists(opath):
        os.makedirs(opath)
    args += ' --outfolder {} '.format(opath)
    odtype = np.dtype([('name', np.unicode, 32), ('ra', np.float32), ('dec', np.float32)])
    if src_arr is not None:
        if isinstance(src_arr, dict):
            src_arr=dict_to_nparray(src_arr, dtype=odtype)
        src_arr = np.atleast_1d(src_arr)

        mask = np.array(['4FGL' not in sname for sname in src_arr['name']])
        if len(src_arr[mask])>0:
            xml_path = os.path.join(opath, 'add_source.xml')
            generate_src_xml(src_arr[mask], xml_path)
            args += '--xml_path {}'.format(xml_path)
    with open('./slurm_draft.xml', 'r') as f:
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

#def get_lc_time(src_of_interest):
#    catalog = fits.open('./gll_psc_v16.fit')
#    nph0 = 0.287966337957
#    names = catalog[1].data['Source_Name']
#    ind = np.where(names==src_of_interest)[0][0]
#    src = catalog[1].data[ind]
#    spec = get_spectrum(src)
#    nph = scipy.integrate.quad(lambda E: spec(E)*E, 1e3, 1e5)[0]*1e4
#    print nph
#    return  28. * (nph0/nph)**(0.6)


def get_lc_time3fgl(src_of_interest, emin=1e3):
    data = fits.open('gll_psc_v16.fit')
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
    data = fits.open('gll_psc_v19.fit')
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
    t_sens = (9/ts_sum)*(365*8)
    print t_sens
    if t_sens < 100:
        t_disc = (25/ts_sum)*(365*8)
        return np.max([t_disc, 56])
    return np.min([t_sens, 200])

#Deprecated
#def get_spectrum(src):
#    '''
#    src is a row from the 3FGL catalog
#    '''
#    print('Create Spectrum for {}'.format(src['Source_Name']))
#    if src['SpectrumType'] == 'PowerLaw':
#        print('{}*(E/{})^-{}'.format(src['Flux_Density'],
#                                     src['Pivot_Energy'],
#                                     src['Spectral_Index']))
#        return lambda x: src['Flux_Density']*(x/src['Pivot_Energy'])**(-src['Spectral_Index'])
#    elif src['SpectrumType'] == 'LogParabola':
#        eq_str='{N}*(E/{piv:.2f})^-({a:.2f}+{b:.2f}*log(E/{piv:.2f})'
#        print eq_str.format(N=src['Flux_Density'], piv=src['Pivot_Energy'],
#                            a=src['Spectral_Index'],b=src['beta'])
#        return lambda x: src['Flux_Density']*(x/src['Pivot_Energy'])**(-(src['Spectral_Index']+src['beta']*np.log(x/src['Pivot_Energy'])))
#    else:
#        print('Unkown Spectrum Type: Fall back to Powerlaw')
#        return lambda x: src['Flux_Density']*(x/src['Pivot_Energy'])**(-src['Spectral_Index'])
#

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
        os.system('convert -density 150 {} -quality 90 {}'.format(filename, filename.replace('.pdf', '.png')))
        images.append(imageio.imread(filename.replace('.pdf', '.png')))
    if len(images)>0:
        imageio.mimsave(opath, images)
    return



def get_add_info(src_of_interest):
    data = fits.open('/scratch9/tglauch/realtime_service/main/gll_psc_v19.fit')
    eflux = data[1].data['Energy_Flux100']
    inds = np.argsort(eflux)[::-1]
    eflux = eflux[inds]
    src_names = data[1].data['Source_Name'][inds]
    if len(np.where(src_names==src_of_interest)[0]) == 0:
        return ''
    src_ind = np.where(src_names==src_of_interest)[0][0]
    ostr = ' | Energy flux (E $\geq$ 100MeV): {:.2e} erg/cm$^2$/s [Top {:.1f}\% in 4FGL]'
    return ostr.format(eflux[src_ind], 1.*src_ind/len(eflux)*100)

def get_sources(ra, dec):
    src_dict = {}
    for key in fields:
        src_dict[key] = []
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
        for i in range(len(files[key]['keys'])):
            app =np.atleast_1d(temp[files[key]['keys'][i]])
            src_dict[fields[i]].extend(app)
    print src_dict
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
            if diff < 0.2:
                if not src_dict['name'][i + j + 1] in src_dict['alt_name'][i] \
                    and src_dict['name'][i + j + 1] != src_dict['name'][i]:
                        src_dict['alt_name'][i].append(src_dict['name'][i + j + 1])
                to_delete.append(i + j + 1)
        for k, j in enumerate(to_delete):
            del(src_dict['name'][j - k])
            del(src_dict['ra'][j - k])
            del(src_dict['dec'][j - k])
        i += 1

    fmt_str = '*{}* {}\n ra: {:.2f} deg |  dec: {:.2f} deg | distance: {:.2f} deg [ra: {:.2f} deg , dec: {:.2f} deg]'
    gcd = [GreatCircleDistance(np.radians(src_dict['ra'][i]), np.radians(src_dict['dec'][i]), np.radians(ra),
           np.radians(dec)) for i in range(len(src_dict['name']))]
    src_dict['dist'] = np.degrees(gcd)
    inds = np.argsort(gcd)
    out_str = '' 
    for i in inds:
        add_info = get_add_info(src_dict['name'][i])
        ostr= fmt_str.format(src_dict['name'][i],
                             add_info,
                             src_dict['ra'][i],
                             src_dict['dec'][i],
                             np.degrees(gcd[i]),
                             src_dict['ra'][i] - ra,
                             src_dict['dec'][i] - dec)
        if len(src_dict['alt_name'][i]) > 0:
            ostr +=  '\n Associations: {}'.format(', '.join(src_dict['alt_name'][i]))
        print(ostr)
        out_str += ostr + '\n\n'
    return dict_to_nparray(src_dict, dtype=odtype), out_str
