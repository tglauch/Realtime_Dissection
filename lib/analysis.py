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
from astropy.io import fits
import warnings
import copy

marker_colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']

src_encoding = {'5BZQ' : -2, '5BZB': -2 , '5BZU': -2, '3HSP': -1, 'CRATES': -3,
                '3FGL': -10, '4FGL': -10, 'FHL': -10}

sc_file_path = '/scratch8/tglauch/spacecraft_files/current.fits'

files = collections.OrderedDict([
     ('4fgl', {'file': '4fgl.1.csv',
              'keys': ['Source_Name', 'RA', 'Dec']}),
     ('3fgl', {'file': '3fgl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('3fhl', {'file': '3fhl.1.csv',
              'keys': ['name', 'ra', 'dec']}),
     ('3hsp', {'file': '3hsp.1.csv',
              'keys': ['Name', 'ra', 'dec']}),
     ('5bzcat', {'file': '5bzcat.1.csv',
                'keys': ['Name', 'RAJ2000', 'DEJ2000']}),
     ('FermiGRB', {'file': 'fgrb.1.csv',
                   'keys' : ['gcn_name','ra','dec']}),
     ('fermi8yr', {'file': 'fermi8yr.1.csv',
                  'keys': ['Source_Name', 'RAJ2000', 'DEJ2000']}),])
#     ('crates', {'file': 'crates.1.csv',
#                'keys': ['name', 'ra', 'dec']})])


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
        self.radius = 200
        self.max_dist = 3.
        self.this_path = None
        return


    def calc_gal_coords(self):
        c = SkyCoord(self.ra, self.dec, frame='icrs', unit="deg")
        self.gal_b = c.galactic.b.deg
        self.gal_l = c.galactic.l.deg
        return


    def make_ts_map_plots(self):
        for tsm in self.ts_maps:
            print('Make plot for {}'.format(tsm))
            try:
                plot.make_ts_plot_legend(tsm, self.srcs, self.max_dist)
                plot.make_ts_plot(tsm, self.srcs, self.vou_sources,
                                  plt_mode='tsmap', error90 = self.err90)
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
            self.radius = np.max([1.5 * 60. * self.err90, 60.])
        elif isinstance(self.err90, Ellipse):
            self.radius = np.max([1.5 * self.err90.get_max_extension(), 60])
        self.max_dist = self.radius/60.
        return
        

    def mask_sources(self):
        if isinstance(self.err90, float):
            for src in self.srcs:
                src.in_err=True
                if src.dist > self.max_dist:
                    src.in_err = False
        elif isinstance(self.err90, Ellipse):
            for src in self.srcs:
                src.in_err = True
                ra_bool = np.abs(src.ra - self.err90.center_ra) > (1.5 * self.err90.ra_ax)
                dec_bool = np.abs(src.dec - self.err90.center_dec) > (1.5 * self.err90.dec_ax)
                if (ra_bool | dec_bool):
                    src.in_err = False
        return

    def make_ts_map(self, ts_map_path, emin=None, trange=None):
        sargs = ' --free_radius {} --data_path {} --use_4FGL --emin {} --ra {} --dec {} --roiwidth {} --free_diff'
        if emin is None:
            emin = self.ts_emin
        if trange is not None:
           sargs = sargs + ' --time_range {} {}'.format(trange[0], trange[1])
        sargs = sargs.format(get_68_psf(emin), self.fermi_data, emin, self.ra, self.dec, 2*self.radius/60.)
        self.ts_maps.append(ts_map_path)
        submit_fit(sargs, ts_map_path, sub_file=self.id+'.sub', ana_type='TS_Map', partition='xtralong')
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
        plot.make_counterparts_plot(self.bpath, self.ra, self.dec, self.radius/60., save_path = self.vou_out,  vou_cand=self.vou_sources, srcs=self.srcs,
                                    max_dist=self.max_dist, legend=False, yaxis=True, error90=self.err90)
        return


    def make_aladin_plot(self):
        with open('./htmlcode/aladin.html') as ifile:
             aladin_code = ifile.read()

        temp_str  = 'overlay.add(A.circle( {ra} ,  {dec}, .05, {{color: \'{color}\'}} )); \n \
                     overlay.add(A.circle( {ra} ,  {dec}, .0027, {{color: \'{color}\'}} )); \n \
                     cat.addSources([A.source({ra} , {dec}, {{name: \' Nr. {name}\'}})]);\n'
        circle_str = ''
        for j, src in enumerate(self.vou_sources):
            circle_str += temp_str.format(ra=src[0], dec=src[1], color='blue', name=j)
        for src in self.srcs:
            circle_str += temp_str.format(ra=src.ra, dec=src.dec, color='red', name=src.name)
        cpath = os.path.join(self.bpath, 'contour.txt')
        cpath_df= os.path.join(self.bpath, 'contour_df.txt')
        nlist = []
        if os.path.exists(cpath):
            cdata = np.degrees(np.genfromtxt(cpath, delimiter=' '))
            for i in cdata:
                nlist.append(list(i)) 
        elif os.path.exists(cpath_df):
            cdata =  np.degrees(np.genfromtxt(cpath_df, delimiter=' '))
            for i in cdata:
                nlist.append([i[1], i[0]])
        elif isinstance(self.err90, Ellipse):
            t = np.linspace(0, 2*np.pi, 100)
            ra = self.err90.center_ra + self.err90.ra_ax * np.cos(t)
            dec = self.err90.center_dec + self.err90.dec_ax * np.sin(t)
            for i in (zip(ra,dec)):
                nlist.append(list(i))
        circle_str += 'overlay.add(A.polyline({pos_arr}, {{color:\'white\'}}));'.format(pos_arr=str(nlist))           
        aladin_code = aladin_code.format(ra=self.ra, dec=self.dec, circles=circle_str)
        return aladin_code
            


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
        for i in range(1): # With new version only need to run ones?!?!
            subprocess.call(cmd)

        # Setup Variables
        out_str = self.get_sources()
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
            self.mode = 'end'
            self.from_gcn(self.gcn, update=True)
        return

    def from_gcn(self, url, update=False):
        gcn_dict = read_gcn(url)
        if gcn_dict['NOTICE_DATE'] != self.notice_date:
            self.ra = gcn_dict['SRC_RA']
            self.dec = gcn_dict['SRC_DEC']
            self.err90 = gcn_dict['SRC_ERROR']
            if 'SRC_ERROR50' in gcn_dict.keys():
                self.err50 = gcn_dict['SRC_ERROR50']
            self.mjd = gcn_dict['MJD']
            #if self.notice_date != 'None':
            #    print('Update GCN and VOU ROI')
            #    self.ROI_analysis()
            self.notice_date = gcn_dict['NOTICE_DATE']
            self.gcn = url
        if update is False:
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
        return out_str

    def set_src_plot_style(self):
        # srcs arleady need an attribute plt_code that is compatible with the cat dict 
        marker_styles = np.unique([src_encoding[key] for key in src_encoding.keys()])
        marker_counter = dict()
        for sty in marker_styles:
            marker_counter[sty] = 0
        for src in self.srcs:
            if not src.in_err:
                src.plt_style = copy.copy(cat_dict[src.plt_code][1])
                src.plt_style['color'] = 'grey'
                continue
            if marker_counter[src.plt_code] > len(marker_colors):
                print('no more colors for this sourcetype {}.Falling back to default'.format(src.name))
                src.plt_style = copy.copy(cat_dict[src.plt_code][1])
                continue
            src.plt_style = copy.copy(cat_dict[src.plt_code][1])
            src.plt_style['color'] = marker_colors[marker_counter[src.plt_code]]
            marker_counter[src.plt_code] += 1
        for src in self.srcs:
            print('name: {} dist: {:.1f} in_err: {} plt_style: {}'.format(src.name, src.dist, src.in_err,
                                                                          src.plt_style))
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
        mask = [src.in_err for src in self.srcs]
        counterparts = np.array([src for src in np.array(self.srcs)[mask] if not 'CRATES' in src.name])
        #inds = np.argsort([src.dist for src in counterparts])
       # counterparts = counterparts[inds]
        html_files = os.path.join(self.bpath, '{}/IceCube{}_files'.format(self.event_name, self.event_name[2:]))
        if os.path.exists(html_files):
            shutil.rmtree(html_files)
        shutil.copytree('./htmlcode/', html_files)
        src_string = ''
        for i, src in enumerate(counterparts):
            with open(os.path.join(html_files, 'src.html'), 'r') as f:
                in_str = f.read()
                src_string += in_str.format(num=i+1,
                                            ra = src.ra,
                                            dec = src.dec,
                                            src_name=src.name,
                                            dist=src.dist,
                                            alt_names= ', '.join(src.names),
                                            sed_path='CandSED{}.png'.format(i+1),
                                            lc_path='CandLC{}.png'.format(i+1))

        c = SkyCoord(self.ra, self.dec, frame='icrs', unit="deg")
        gal = c.galactic      
        aladin_str = self.make_aladin_plot()
        with open(os.path.join(html_files, 'template.html'), 'r') as f:
            code = f.read()
            code = code.format(mjd = self.mjd,
                               src_summary=src_string,
                               ev_name_full = self.event_name.replace('IC', 'IceCube-')+'A',
                               ev_name = self.event_name,
                               dec = self.dec,
                               ra = self.ra,
                               gal_lat = gal.b.deg,  
                               ts_mjd1 = self.tsmjd1,
                               ts_mjd2 =  self.tsmjd2,
                               ts_full_mjd1 = self.mjd_range[0],
                               ts_full_mjd2 = self.mjd_range[1],
                               ts_length = self.tsmjd2 - self.tsmjd1,
                               aladin_text = aladin_str,
                               )
            
        with open(os.path.join(html_files, 'index.html'), 'w') as f:
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
        if os.path.exists(os.path.join(self.bpath, 'ts_map/legend.png')):
            shutil.copyfile(os.path.join(self.bpath, 'ts_map/legend.png'),
                            os.path.join(html_files, 'legend.png'))
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
            print self.mode
            if len(folders) >0:
                if self.mode == 'end':
                    sed_path = os.path.join(lc_base,sorted(folders)[-1])
                else:
                    for fold in folders:
                        if (self.mjd <= float(fold.split('_')[1])) and (self.mjd > float(fold.split('_')[0])):
                            sed_path = os.path.join(lc_base,fold)
                            break
                print(os.path.join(html_files, 'CandSED{}.png'.format(i+1)))
                if os.path.exists(os.path.join(sed_path, 'sed.png')):
                    shutil.copyfile(os.path.join(sed_path, 'sed.png'), os.path.join(html_files, 'CandSED{}.png'.format(i+1)))
            else:
                print('SED seems to be not ready, yet')
            if os.path.exists(lc_path):
                shutil.copyfile(lc_path, os.path.join(html_files, 'CandLC{}.png'.format(i+1)))
        self.upload_html(html_files)
        return

    def upload_html(self, html_files):
        cmd = 'rsync -r {}/* tglauch@cobalt.icecube.wisc.edu:/home/tglauch/public_html/Reports/{}/'.format(html_files, self.event_name)
        print(cmd)
        os.system(cmd)
        return
