from functions import GreatCircleDistance, path_settings, vou_path, get_68_psf, submit_fit, cat_dict
import numpy as np
from source_class import Source, Ellipse 
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
import plot
import shutil
import pyfits as fits
import warnings
import copy

marker_colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']

src_encoding = {'5BZQ' : -2, '5BZB': -2 , '5BZU': -2, '3HSP': -1, 'CRATES': -3,
                '3FGL': -10, '4FGL': -10, 'FHL': -10}

sc_file_path = '/scratch9/tglauch/Realtime_Dissection/sc_files/current.fits'

files = collections.OrderedDict([
     ('4fgl', {'file': '4fgl.1.csv',
              'keys': ['Source_Name', 'RA', 'Dec']}),
     ('3fgl', {'file': '3fgl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('3fhl', {'file': '3fhl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('FermiGRB', {'file': 'fgrb.1.csv',
                   'keys' : ['GRB','RAJ2000','DEJ2000']}),
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
        if self.mjd is not None:
            t = Time.now()
            if np.abs(self.mjd - t.mjd)<2:
                self.mode = 'end'
            else:
                self.mode = 'mid'
        self.mjd_range = None
        self.gcn = 'No GCN given'
        self.vou_out = None
        self.fermi_data = None
        self.ts_maps = []
        self.srcprob_path = None
        self.id = id_generator(size=5) 
        self.notice_date = 'None'
        self.radius = 180
        self.this_path = None
        return


    def calc_gal_coords(self):
        c = SkyCoord(self.ra, self.dec, frame='icrs', unit="deg")
        self.gal_b = c.galactic.b.deg
        self.gal_l = c.galactic.l.deg
        return


    def make_ts_map_plots(self):
        for tsm in self.ts_maps:
            try:
                plot.make_ts_plot(tsm, self.srcs, self.vou_sources,
                                  plt_mode='tsmap', error90 = self.err90)
                plot.make_ts_plot_legend(tsm, self.srcs, self.max_dist)
            except Exception as inst:
                warnings.warn("Couldn't create TS map {}".format(tsm))
                print(inst)
            try:
                plot.make_ts_plot(tsm, self.srcs, self.vou_sources,
                                  plt_mode='residmap',  error90=self.err90)
            except Exception as inst:
                warnings.warn("Couldn't create residual map {}".format(tsm))
                print(inst)
        return

    def calc_src_probs(self, srcprob_path, emin=None):
        sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --free_diff'
        if emin is None:
            emin = self.ts_emin
        sargs = sargs.format(get_68_psf(emin), self.fermi_data, emin, self.ra, self.dec)
        self.srcprob_path = srcprob_path
        submit_fit(sargs, srcprob_path, srcs=self.srcs, sub_file=self.id+'.sub', ana_type='srcprob', partition='xtralong')
        return

    def adaptive_radius(self):
        if isinstance(self.err90, float):
            self.radius = np.max([2 * 60. * self.err90, 60.])
        elif isinstance(self.err90, Ellipse):
            self.radius = np.max([2 * self.err90.get_max_extension(), 60])
        return
        

    def make_ts_map(self, ts_map_path, emin=None, trange=None):
        sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --roiwidth {} --free_diff'
        if emin is None:
            emin = self.ts_emin
        if trange is not None:
           sargs = sargs + ' --time_range {} {}'.format(trange[0], trange[1])
        sargs = sargs.format(get_68_psf(emin), self.fermi_data, emin, self.ra, self.dec, 2*self.radius/60.)
        self.ts_maps.append(ts_map_path)
        submit_fit(sargs, ts_map_path, sub_file=self.id+'.sub', ana_type='TS_Map', partition='xtralong') # change back to xtralong
        return

    def make_pdf(self, src_latex, final_pdf=True):
        with open(os.path.join(self.bpath, 'vou_blazar/full_output'), 'r') as f:
            full_out = f.read().encode('utf8').\
                replace('\x1b', '').replace('nu_p', 'nu peak')
            full_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', full_out)
            full_out = re.sub('(Match[^\\\\]*)', r'\\textbf{\1}', full_out)
            full_out = re.sub('(Dist[^\\\\]*)', r'\\textbf{\1}', full_out)
            full_out = [i.strip() + ' \n' for i in full_out.split('\n')]
            full_out = [i if i != '\\\\ \n' else '\\\\\\\ \n' for i in full_out]
            full_out = ' '.join(full_out)
        with open(os.path.join(self.bpath, 'vou_blazar/short_output'), 'r') as f:
            short_out = f.read().encode('utf8').replace('\x1b', '')
            short_out = re.sub('\*([^\*]*)\*', r'\\textbf{\1}', short_out)
        with open(os.path.join(self.this_path, 'latex','template.tex')) as f:
            template = f.read()

        c = SkyCoord(self.ra, self.dec, frame='icrs', unit="deg")
        gal = c.galactic
        if not final_pdf:
            prelim = '{\color{red} Warning: The following results are all preliminary and will be continously updated whenever calculations are finished.}'
        else:
            prelim = ' '
        t = Time.now()
        t_now_str = t.iso
        out = template.format(prelim=prelim, date=self.mjd, ra=self.ra, dec=self.dec,
                              emin=1.*self.emin/1000., l=gal.l.deg, b=gal.b.deg,
                              cat_srcs=short_out,
                              rx_map=os.path.join(self.vou_out, 'RX_map.eps'),
                              vou_pic=os.path.join(self.vou_out, 'counterparts.png'),
                              ts_map=os.path.join(self.bpath, 'ts_map/tsmap.png'),
                              ts_map_short=os.path.join(self.bpath, 'ts_map_short/tsmap.png'),
                              ts_map_legend=os.path.join(self.bpath, 'ts_map_short/legend.png'),
                              vou_output=full_out, event=self.event_name,
                              tsmjd1=self.tsmjd1, tsmjd2=self.tsmjd2, src_latex=src_latex,
                              mjd1=self.mjd_range[0], mjd2=self.mjd_range[1],
                              energy=self.emin/1000., tsemin = self.ts_emin/1e3,
                              gcnurl = self.gcn, createdon=t_now_str)

        latex_path = os.path.join(self.bpath, self.event_name + '.tex')
        if os.path.exists(latex_path):
            os.remove(latex_path)
        with open(latex_path, 'w+') as f:
            f.write(out)
        if not os.path.exists(os.path.join(self.bpath, 'sample.bib')):
            shutil.copyfile('./latex/sample.bib', os.path.join(self.bpath, 'sample.bib'))
        os.chdir(self.bpath)
        cmds  = [
            ['pdflatex', '-interaction', 'nonstopmode', self.event_name + '.tex'],
            ['bibtex', self.event_name + '.aux'],
            ['pdflatex', '-interaction', 'nonstopmode',  self.event_name + '.tex'],
            ['pdflatex', '-interaction', 'nonstopmode',  self.event_name + '.tex'],
        ]
        for c in cmds:
            subprocess.call(c)
        os.chdir(self.this_path)
        self.pdf_out_path = os.path.join(self.bpath, self.event_name + '.pdf') 
        return


    def make_counterparts_plot(self):
        self.get_vou_candidates()
        plot.make_counterparts_plot(self.ra, self.dec, save_path = self.vou_out,  vou_cand=self.vou_sources, srcs=self.srcs,
                                    max_dist=self.max_dist, legend=False, yaxis=True, error90=self.err90)
        return


    def ROI_analysis(self):
        self.vou_out = os.path.join(self.bpath, 'vou_blazar')
        if not os.path.exists(self.vou_out):
            os.makedirs(self.vou_out)
        os.chdir(self.vou_out)

        if isinstance(self.err90, float):
            print('Run VOU with cicular error')
            cmd = [vou_path, str(self.ra), str(self.dec), str(self.radius), str(self.err90 * 60)]
        elif isinstance(self.err90, Ellipse):
            print('Run VOU with elliptical error')
            cmd = [vou_path]
            cmd.extend(self.err90.get_vou_cmd(self.radius))
        else:
            print('This type of error regions is not supported yet')
        print cmd
        for i in range(2):
            subprocess.call(cmd)

        # Setup Variables
        out_str = self.get_sources()
        self.make_counterparts_plot()
        self.calc_gal_coords()
        headline = '*Result for {}* [ra: {:.1f}  dec: {:.1f} , gal_l: {:.1f}, gal_b: {:.1f} ]\n'.format(self.event_name,
                                                                               self.ra, self.dec, self.gal_l, self.gal_b)
        sum_text = headline + out_str
        print_to_slack(sum_text)

        # Convert VOU Blazar Output
        rx_ps = os.path.join(self.vou_out, 'RX_map.ps')
        os.system('ps2eps -B ' + rx_ps)
        cand_ps = os.path.join(self.vou_out, 'candidates.ps')
        os.system('ps2eps -B ' + cand_ps)
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
        os.chdir(self.this_path)
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
        MET = get_data(self.ra, self.dec, sc_file = sc_file_path , **args)
        MJD = [MET_to_MJD(float(i)) for i in MET]
        self.mjd_range = MJD
        return

    def update_gcn(self):
        if self.gcn != 'No GCN given':
            self.from_gcn(self.gcn)
        return

    def from_gcn(self, url):
        gcn_dict = read_gcn(url)
        if gcn_dict['NOTICE_DATE'] != self.notice_date:
            self.ra = gcn_dict['SRC_RA']
            self.dec = gcn_dict['SRC_DEC']
            self.err90 = gcn_dict['SRC_ERROR']
            if 'SRC_ERROR50' in gcn_dict.keys():
                self.err50 = gcn_dict['SRC_ERROR50']
            self.mjd = gcn_dict['MJD']
            if self.notice_date != 'None':
                print('Update GCN and VOU ROI')
                self.ROI_analysis()
            self.notice_date = gcn_dict['NOTICE_DATE']
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


    def sort_sources_by_prio(self, prio=['3HSP', '5BZB', '5BZQ',' 5BZU', '3FGL', '4FGL', 'FHL', 'CRATES']):
        new_order = []
        print self.srcs
        for sclass in prio:
            for i in range(len(self.srcs)):
                if i in new_order:
                    continue
                if (sclass in self.srcs[i].name) or (np.any([True if sclass in j else False for j in
                    self.srcs[i].names])):
                        new_order.append(i)
                        self.srcs[i].plt_code = src_encoding[sclass] 
        new_order = np.array(new_order)
        self.srcs = [self.srcs[j] for j in new_order]
        # Remove Later!!!
        self.set_src_plot_style()
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
                if key == 'FermiGRB':
                    name = 'GRB' + src[files[key]['keys'][0]]
                else:
                    name = src[files[key]['keys'][0]]
                src_obj = Source(name,
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
        self.sort_sources_by_prio()
        out_str = ''
        for src in self.srcs:
            add_info = self.get_add_info(src.name)
            ostr= fmt_str.format(src.name, add_info, src.ra, src.dec,
                                 src.dist, src.ra - self.ra, src.dec - self.dec)
            if len(src.names) > 0:
                ostr +=  '\n Associations: {}'.format(', '.join(src.names))
            print(ostr)
            out_str += ostr + '\n\n'
        self.set_src_plot_style()
        return out_str

    def set_src_plot_style(self):
        # srcs arleady need an attribute plt_code that is compatible with the cat dict 
        marker_styles = np.unique([src_encoding[key] for key in src_encoding.keys()])
        marker_counter = dict()
        for sty in marker_styles:
            marker_counter[sty] = 0
        for src in self.srcs:
            if src.dist > self.max_dist:
                src.plt_style = copy.copy(cat_dict[src.plt_code][1])
                continue
            if marker_counter[src.plt_code] > len(marker_colors):
                print('no more colors for this sourcetype {}.Falling back to default'.format(src.name))
                src.plt_style = copy.copy(cat_dict[src.plt_code][1])
                continue
            src.plt_style = copy.copy(cat_dict[src.plt_code][1])
            src.plt_style['color'] = marker_colors[marker_counter[src.plt_code]]
            marker_counter[src.plt_code] += 1
            print src.plt_style
        return    


    def get_vou_candidates(self):
        x = np.genfromtxt(os.path.join(self.vou_out, 'find_out_temp.txt'),
                          dtype=[np.float, np.float, int])
        mask = []
        for i in x:
            if i[2] < 0 :
                 mask.append(False)
            else:
                mask.append(True)
        self.vou_sources = x[mask]
        return 
   
    def get_add_info(self, src_of_interest):
        data = fits.open(os.path.join(self.this_path,'lib',  'gll_psc_v19.fit'))
        eflux = data[1].data['Energy_Flux100']
        inds = np.argsort(eflux)[::-1]
        eflux = eflux[inds]
        src_names = data[1].data['Source_Name'][inds]
        if len(np.where(src_names==src_of_interest)[0]) == 0:
            return ''
        src_ind = np.where(src_names==src_of_interest)[0][0]
        ostr = ' | Energy flux (E $\geq$ 100MeV): {:.2e} erg/cm$^2$/s [Top {:.1f}\% in 4FGL]'
        return ostr.format(eflux[src_ind], 1.*src_ind/len(eflux)*100)


    def create_html(self):
        mask = [True if src.dist < self.max_dist else False for src in self.srcs]
        counterparts = np.array([src for src in np.array(self.srcs)[mask] if not 'CRATES' in src.name])
        inds = np.argsort([src.dist for src in counterparts])
        counterparts = counterparts[inds][:2]
        html_files = os.path.join(self.bpath, '{}/IceCube{}_files'.format(self.event_name, self.event_name[2:]))
        if os.path.exists(html_files):
            shutil.rmtree(html_files)
        shutil.copytree('./htmlcode/IceCubeTemplate_Candidates_files/', html_files)
        index_html = os.path.join(self.bpath, '{}/index.html'.format(self.event_name))
        shutil.copyfile('./htmlcode/IceCubeTemplate_{}Candidates.html'.format(len(counterparts)),
                        index_html)
        with open(index_html, 'r') as f:
            code = f.read()
            code = code.replace('IceCube yyyy', 'IceCube {}'.format(self.event_name[2:]))
            code = code.replace('ICyyyy', self.event_name)
            code = code.replace('MJD aaaa', 'MJD {}'.format(self.mjd))
            code = code.replace('gggg', str(9999))
            code = code.replace('tttt', self.event_name[2:])
            code = code.replace('rrrr', self.event_name[2:])
            
        with open(index_html, 'w') as f:
            f.write(code)
        
        rxmap_eps = os.path.join(self.vou_out, 'RX_map.eps')
        rxmap_png = os.path.join(self.vou_out, 'RX_map.png')
        os.system('convert -density 288 {} {}'.format(rxmap_eps, rxmap_png))

        #candidates_eps = os.path.join(self.vou_out, 'candidates.eps')
        candidates_png = os.path.join(self.vou_out, 'counterparts.png')
        #os.system('convert -density 288 {} {}'.format(candidates_eps, candidates_png))
        if os.path.exists(rxmap_png):
            shutil.copyfile(rxmap_png,
                            os.path.join(html_files, 'VOU-RXmap.png'))
        if os.path.exists(candidates_png):
            shutil.copyfile(candidates_png,
                            os.path.join(html_files, 'VOU-candidates.png'))
        if os.path.exists(os.path.join(self.bpath, 'ts_map/tsmap.png')):
            shutil.copyfile(os.path.join(self.bpath, 'ts_map/tsmap.png'), 
                            os.path.join(html_files, 'tsmap-full.png'))
        if os.path.exists(os.path.join(self.bpath, 'ts_map_short/tsmap.png')):
            shutil.copyfile(os.path.join(self.bpath, 'ts_map_short/tsmap.png'),
                            os.path.join(html_files, 'tsmap_200.png'))
        pdf_path = os.path.join(self.bpath, self.event_name + '.pdf')
        if os.path.exists(pdf_path):
            shutil.copyfile(pdf_path, os.path.join(html_files, self.event_name + '.pdf'))
        for i, src in enumerate(counterparts):
            print i, src.name
            lc_base = src.lc_path
            print lc_base
            lc_path = os.path.join(lc_base, 'lightcurve.png')
            folders = [fold for fold in os.listdir(lc_base) if os.path.isdir(os.path.join(lc_base, fold))]
            if len(folders) >  0:
                if self.mode == 'end':
                    sed_path = os.path.join(lc_base,sorted(folders)[-1])
                else:
                    for fold in folders:
                        if (self.mjd <= float(fold.split('_')[1])) and (self.mjd > float(fold.split('_')[0])):
                            sed_path = os.path.join(lc_base,fold)
                            break
                print(os.path.join(html_files, 'CANDSED{}.png'.format(i+1)))
                shutil.copyfile(os.path.join(sed_path, 'sed.png'), os.path.join(html_files, 'CandSED{}.png'.format(i+1)))
            else:
                print('SED seems to be not ready, yet')
            print('CANDLC{}.png'.format(i+1))
            if os.path.exists(lc_path):
                shutil.copyfile(lc_path, os.path.join(html_files, 'CandLC{}.png'.format(i+1)))
        self.upload_html()
        return

    def upload_html(self):
        os.system('scp -r {} tglauch@cobalt.icecube.wisc.edu:/home/tglauch/public_html/Reports/'.format(os.path.join(self.bpath, self.event_name)))
        return
