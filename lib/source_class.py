import numpy as np
from functions import GreatCircleDistance, vou_path, path_settings, get_lc_time3fgl, \
                      get_lc_time4fgl,get_68_psf, submit_fit, make_gif
from astropy.io import fits
from myfunctions import MET_to_MJD, dict_to_nparray, ts_to_pval, pval_to_sigma
import os
import plot
import warnings
import shutil
from add_classes import Lightcurve, Ellipse
import multiprocessing
from read_catalog import read_from_observation
import requests
 
ul_ts_threshold = 4
MeV_to_erg = 1.60218e-6

class Source(object):
    """Class holding all the information for a Counterpart Candidate"""

    def __init__(self, name, ra, dec):
        self.name = name
        self.names = []
        self.ra = ra  
        self.dec = dec
        self.lightcurves = dict()
        self.seds = dict()
        self.mw_idata = None
        self.ovro = None
        self.swift = None

    def make_sed_lightcurve(self, lcs=['default']):
        print('Make SED lightcurve for {}'.format(self.name))
        print self.lightcurves.keys()
        main_lc = self.lightcurves[lcs[0]]
        #if not hasattr(self, 'mw_idata'):
        self.set_mw_data()
        pool = multiprocessing.Pool()
        y_vals = []
        if self.mw_idata is not None:
            inds = (self.mw_idata[:,1] > 0) & (self.mw_idata[:,0] < 1e22)
            y_vals.extend(self.mw_idata[:,1][inds])

        for basepath in main_lc.time_window_results:
            if os.path.exists(os.path.join(basepath, 'sed.npy')):
                sed = np.load(os.path.join(basepath, 'sed.npy'), allow_pickle=True)[()]
            else:
                continue
            m = sed['ts'] < ul_ts_threshold
            y_vals.extend(sed['e2dnde'][~m]*MeV_to_erg)
            y_vals.extend(sed['e2dnde_ul95'][m]*MeV_to_erg)

        if len(y_vals) ==0:
            y_min = 1e-15
            y_max = 1e-9
        else:
            y_max = 2 * np.min([1e-8, 1.1 * np.max(y_vals)])
            y_min = 0.5 * np.min([0.9 * np.min(y_vals), 1e-15])
       
        for i in range(len(main_lc.time_windows)):
            try:
                seds_list = [(main_lc.time_window_results[i], 'k', 'red' , True, True, True),
                             (self.seds['default'], 'grey','grey', True, True, True)]
            #    for j in range(1,len(lcs)):    
            #        seds_list.append((self.lightcurves[lcs[j]].time_window_results[i],
            #                          'k', 'blue' , True, True, False))
                kwarg_dict = {'mw_idata': self.mw_idata, 'dec':self.dec, 'twindow':main_lc.time_windows[i], 'y_min':y_min,
                              'y_max': y_max, 'add_text':True}
                pool.apply_async(plot.make_sed_plot, (seds_list,), kwarg_dict) 
            except Exception as inst:
                warnings.warn("Couldn't create SED for source {}".format(self.name))
                print(inst)
        pool.close()
        pool.join()

        #Create all year SED
        try:
            plot.make_sed_plot([(self.seds['default'], 'grey', 'grey', True, True, True)],
                               mw_idata=self.mw_idata, dec=self.dec, y_min=y_min, y_max=y_max, add_text=True)
        except Exception as inst:
                warnings.warn("Couldn't create all year SED")
                print(inst)

        # Create GIF
        try:
            make_gif(main_lc.bpath)
        except Exception as inst:
            warnings.warn("Couldn't create an animated light curve {}".format(self.name))
            print(inst)
        return 

    def get_ovro_data(self):
        self.ovro = None
        mask = [True if '3FGL' in i else False for i in self.names]
        fgl_name = np.array(self.names)[mask]
        if len(fgl_name)  >0 :
            fgl_name = fgl_name[0]
            print('Looking for OVRO data for {}'.format(self.name))
            ovro_name = fgl_name[5:10] + fgl_name[12:]
            html= requests.get('http://www.astro.caltech.edu/ovroblazars/HL28/csv/{}.csv'.format(ovro_name))
            if html.ok:
                print('Data Found')
                with open(os.path.join(self.bpath, 'add_data', 'ovro.csv'), 'w+') as ofile:
                    ofile.write(html.text)
                self.ovro = os.path.join(self.bpath, 'add_data', 'ovro.csv')
            else:
                print('No Data Found!')
        return                 


    def make_lc_plot(self, mjd):
        self.get_ovro_data()
        self.swift=None
        self.collect_xray()
        for i, lc_key in enumerate(self.lightcurves.keys()):
            try:
                if i == 0:
                    plot.make_lc_plot(self.lightcurves[lc_key].bpath, mjd, self.name, radio=self.ovro, xray=self.swift)
                else:
                    plot.make_lc_plot(self.lightcurves[lc_key].bpath, mjd, self.name)
            except Exception as inst:
                warnings.warn("Couldn't create {} lightcurve for {}".format(lc_key, self.name))
                print(inst)
        return 
    
    def angular_distance(self, ra, dec):
        val = GreatCircleDistance(np.radians(self.ra), np.radians(self.dec),
                                  np.radians(ra), np.radians(dec))
        return np.degrees(val)
    
    def angular_distance_to_source(self, src2):
        val = GreatCircleDistance(np.radians(self.ra), np.radians(self.dec),
                                  np.radians(src2.ra), np.radians(src2.dec))
        return np.degrees(val)


    def setup_folders(self, bpath_src):
        self.bpath = bpath_src
        if not os.path.exists(self.bpath):
            os.makedirs(self.bpath)
        self.sed_path = os.path.join(bpath_src, path_settings['sed'])
        if not os.path.exists(self.sed_path):
            os.makedirs(self.sed_path)
        self.lc_path = os.path.join(bpath_src, path_settings['lc'])
        if not os.path.exists(self.lc_path):
            os.makedirs(self.lc_path)
        self.mw_data_path = os.path.join(bpath_src, 'sed.txt')
        if not os.path.exists(os.path.join(bpath_src, 'add_data')):
            os.makedirs(os.path.join(bpath_src, 'add_data'))
        return


    def get_mw_data(self):
        this_path = os.getcwd()
        if '4FGL' in self.name:
            data = fits.open('./lib/gll_psc_v19.fit')[1].data
            ind = np.where(data['Source_Name']==self.name)[0]
            m_ax = float(data[ind]['Conf_95_SemiMajor']) * 60
            min_ax = float(data[ind]['Conf_95_SemiMinor']) * 60
            loc_str= '{} {} {} {}'.format(2 * np.max([m_ax, min_ax]), m_ax, min_ax,
                                          float(data[ind]['Conf_95_PosAng']))
            ofolder = os.path.join(self.bpath, 'vou_counterpart')
            if os.path.exists(ofolder):
                shutil.rmtree(ofolder)
            os.makedirs(ofolder)
            os.chdir(ofolder)
            for i in range(2):
                os.system('{vou_path} {ra} {dec} {loc_str}'.format(vou_path=vou_path, ra=self.ra,
                                                                   dec=self.dec,loc_str=loc_str))
            if not os.path.exists('candidates.eps'):
                fname_new = os.path.join(ofolder, self.name.replace(' ', '_').replace('.','_') + '.ps')
                cand_ps = os.path.join(ofolder, 'candidates.ps')
                os.rename(cand_ps, fname_new)
                os.system('ps2eps -B ' + fname_new)
            else:
                cand_eps = os.path.join(ofolder, 'candidates.eps')
                fname_new = os.path.join(ofolder, self.name.replace(' ', '_').replace('.','_') + '.eps')
                os.rename(cand_eps, fname_new)
            os.chdir(this_path)
            cp_info = np.array(np.genfromtxt(os.path.join(ofolder, 'candidates_posix.txt'), dtype=float), ndmin=2)
            if len(cp_info) == 0:
                cp_ra = self.ra
                cp_dec = self.dec
            else:
                try:
                    cp_ra = cp_info[:,0][0]
                    cp_dec = cp_info[:,1][0]
                    print('Found Counterpart at ra {}, dec {}'.format(cp_ra, cp_dec))
                except Exception as inst:
                    cp_ra = self.ra
                    cp_dec = self.dec
        else:
            cp_ra = self.ra
            cp_dec = self.dec
        pre_str = '{vou_path} {ra} {dec} {loc_str} -s'
        form_str = pre_str.format(vou_path=vou_path, ra=cp_ra, dec=cp_dec,loc_str=1)
        print('Fetch SED with command \n {}'.format(form_str))
        for i in range(2):
            os.system(form_str)
        if os.path.exists('Sed.txt'):
            print('Move Sed.txt to {}'.format(self.mw_data_path))
            shutil.move('Sed.txt', self.mw_data_path)
        else:
            print('There is no File called Sed.txt')
        self.set_mw_data()
        self.get_ovro_data()
        return

    def collect_xray(self):
        if not os.path.exists(self.mw_data_path):
            self.get_mw_data() 
        try:
            idata = np.genfromtxt(self.mw_data_path, skip_header=4,
                                  usecols=[1,2,3,4,6], dtype=[np.float,np.float,np.float,np.float,object])
        except Exception as inst:
            print(inst)
            self.swift=None
            return
        idata = np.array([list(i) for i in idata])
        try:
            idata = idata[idata[:,4] == 'OUSXB']
            self.swift=None
            return
        except Exception as inst:
            print(inst)
            self.swift=None
            return
        idata = np.array(idata[:,0:4],dtype=np.float)
        idata[:,1] = idata[:,1] - idata[:,0]
        idata[:,2] = idata[:,0] - idata[:,2]
        if not os.path.exists(os.path.join(self.bpath, 'add_data')):
            os.makedirs(os.path.join(self.bpath, 'add_data'))
        if len(idata) >0:
            np.savetxt(os.path.join(self.bpath, 'add_data', 'swift.csv'), idata)
            self.swift = os.path.join(self.bpath, 'add_data', 'swift.csv')
        else:
            self.swift = None  
        return


    def set_mw_data(self):
        print('Read multi-wavelength data from {}'.format(self.mw_data_path))
        try:
            self.mw_idata = np.atleast_2d(np.genfromtxt(self.mw_data_path, skip_header=4,
                                          usecols=(0,1,2,3,4)))
        except Exception as inst:
            print(inst)
            self.mw_idata = None
            return
        if len(np.squeeze(self.mw_idata)) == 0:
            self.mw_idata = None
        add_data_path = os.path.join(self.bpath, 'add_data')
        if os.path.exists(add_data_path):
            files = [os.path.join(add_data_path, i) for i in os.listdir(add_data_path)]
            for f in files:
                print('Read from {}'.format(f))
                self.mw_idata = np.concatenate([self.mw_idata, read_from_observation(f)])
        return

    def source_summary(self, bpath, mjd, mode='mid'):
        sed_path = None
        bpath_src = self.bpath
        lc_base = self.lc_path
        lc_path = os.path.join(lc_base, 'lightcurve.pdf')
        folders = [fold for fold in os.listdir(lc_base) if os.path.isdir(os.path.join(lc_base, fold))]
        if len(folders) == 0:
            return ''
        if mode == 'end':
            sed_path = os.path.join(lc_base,sorted(folders)[-1])
        else:
            for fold in folders:
                if (mjd <= float(fold.split('_')[1])) and (mjd > float(fold.split('_')[0])):
                    sed_path = os.path.join(lc_base,fold)
                    break
        if os.path.exists(os.path.join(self.sed_path, 'llh.npy')):
            sed_full_res = np.load(os.path.join(self.sed_path, 'llh.npy'), allow_pickle=True)[()]
            ts = sed_full_res['sources'][self.name]['ts']
            sigma = np.max([0, pval_to_sigma(ts_to_pval(ts, 1))])
            if sigma > 5:
                sigma = np.sqrt(sed_full_res['sources'][self.name]['ts'])
        else:
            sigma = -1
        if sed_path == None:
            return ''
        l_str ='\subsection{{{srcinfo}}}'
        if os.path.exists(os.path.join(sed_path, 'llh.npy')):
            fit_res = np.load(os.path.join(sed_path, 'llh.npy'), allow_pickle=True)[()]
            t_str = '{src_name} $|$ \\small\\textnormal{{{src_info}}}'
            t_str_info = 'ra = {ra:.2f}$^\circ$, dec = {dec:.2f}$^\circ$, $\Sigma$= {sigma:.1f} $\sigma$, $\Delta\psi$ = {dist:.2f}$^\circ$'
            t_str_info = t_str_info.format(sigma=sigma, ra=self.ra, dec=self.dec, dist=self.dist)
            t_str = t_str.format(src_name = self.name, src_info=t_str_info)
            l_str = l_str.format(srcinfo=t_str) # srcname=src['name'], srcinfo=t_str)
        else:
            t_str = '{src_name} $|$ \\small\\textnormal{{{src_info}}}'
            t_str_info = 'ra = {ra:.2f}$^\circ$, dec = {dec:.2f}$^\circ$, $\Delta\psi$ = {dist:.2f}$^\circ$'
            t_str_info = t_str_info.format(ra=self.ra, dec=self.dec, dist=self.dist)
            t_str = t_str.format(src_name = self.name, src_info=t_str_info)
            l_str = l_str.format(srcinfo=t_str)

        l_str += 'Associations: {} \\\\ \n'.format(', '.join(self.names))
        sed_pdf = os.path.join(sed_path, 'sed.pdf')
        try:
            srcprob = fits.open(os.path.join(bpath,'srcprob/ft1_srcprob_00.fits'))
            energy = srcprob[1].data['ENERGY']
            mjd = MET_to_MJD(srcprob[1].data['TIME'])
            src_prob = srcprob[1].data[self.name]
            prob_mask = src_prob > 0.90
            ind = np.argsort(energy[prob_mask])[::-1]
            with open('./latex/tab_template.tex', 'r') as f:
                tab_str = f.read()
            prob_str = ''
            for i in range(np.min([5, len(mjd[prob_mask][ind])])):
                prob_str += '{:.2f} & {:.2f} & {:.2f} \\\\ \n'.format(mjd[prob_mask][ind[i]],
                                                                      src_prob[prob_mask][ind[i]]*100,
                                                                      energy[prob_mask][ind[i]]/1e3)
            tab_str = tab_str.format(events=prob_str)
            l_str += tab_str
        except Exception as inst:
            warnings.warn("Could not find source probabilities")
            print(inst)
        with open('./latex/source_summary.tex', 'r') as infile:
            fig_str = infile.read()
        if os.path.exists(sed_pdf):
            cap_str = 'SED for {}. See the description in section \\ref{{sec:sed}} for more details.'
            l_str += fig_str.format(width = 0.8, path = sed_pdf, caption=cap_str.format(self.name))
        if os.path.exists(lc_path):
            l_str += fig_str.format(width = 0.8, path = lc_path, caption='Light curve for {}'.format(self.name))
        gev_lc = os.path.join(lc_base+'_1GeV', 'lightcurve.pdf')
        if os.path.exists(gev_lc):
            l_str += fig_str.format(width = 0.8, path = gev_lc, caption='1GeV light curve for {}'.format(self.name))
        if ('4FGL' in self.name) | ('3FGL' in self.name):
            cp_candidate_path = os.path.join(bpath_src, 'vou_counterpart', self.name.replace(' ', '_').replace('.','_') + '.eps')
            if os.path.exists(cp_candidate_path):
                l_str += fig_str.format(width = 0.7, path = cp_candidate_path, caption='Possible couterparts for {}'.format(self.name))
        l_str += '\\clearpage \n'
        return l_str


    def make_fixed_binning_lightcurve(self, emin, fermi_data_path, mjd_range, mjd=None, dt_lc=None, mode='end',
                                      name='', add_srcs=None, job_id='fermi'):
        sargs = '--target_src {} --free_radius {} --data_path {} --use_4FGL --emin {} --free_diff'
        sub_args = sargs.format(self.name.replace(' ', '_'), get_68_psf(emin), fermi_data_path, emin)
        if '3FGL' in self.name:
            dt_lc = get_lc_time3fgl(self.name, emin=emin)
        elif '4FGL' in self.name:
            dt_lc = get_lc_time4fgl(self.name, emin=emin)
        else:
            dt_lc = dt_lc

        if mode == 'end':
            time_windows = [[k - dt_lc, k] for k in
                             np.abs(np.arange(-mjd_range[1], -mjd_range[0], dt_lc))]
        elif mode == 'mid':
            time_windows = [[k - dt_lc, k] for k in
                             np.abs(np.arange(mjd+3*dt_lc/2, mjd_range[1], dt_lc))]
            time_windows2 = [[k - dt_lc, k] for k in
                             np.abs(np.arange(-mjd-dt_lc/2, -mjd_range[0], dt_lc))]
            time_windows.extend(time_windows2)

        if (time_windows[-1][1] - mjd_range[0]) < dt_lc:
            del time_windows[-1]
        if name != '':
            ofolder = self.lc_path + '_' + name
        else:
            ofolder = self.lc_path
        this_lc = Lightcurve(ofolder, time_windows)
        this_lc.emin = emin
        this_lc.fermi_data_path = fermi_data_path
        this_lc.mjd = mjd
        this_lc.mjd_range = mjd_range
        this_lc.mode = mode 
        for t_window in time_windows:
            add_str = '{:.1f}_{:.1f}'.format(t_window[0], t_window[1])
            opath = os.path.join(ofolder, add_str)
            this_lc.time_window_results.append(opath)
            submit_fit(sub_args, opath, srcs=add_srcs, sub_file=job_id + '.sub',
                       trange=t_window, partition='kta')
        if name == '':
            self.lightcurves['default'] = this_lc
        else:
            self.lightcurves[name] = this_lc
        return


    def make_sed(self, emin, fermi_data_path, name='', add_srcs=None, job_id='fermi'):
        sargs = '--target_src {} --free_radius {} --data_path {} --use_4FGL --emin {} --free_diff'
        if name != '':
            opath = self.sed_path + '_' + name
        else:
            opath = self.sed_path
        sub_args = sargs.format(self.name.replace(' ', '_'), get_68_psf(emin), fermi_data_path, emin)
        submit_fit(sub_args, opath, srcs=add_srcs, sub_file=job_id + '.sub', trange='', partition='xtralong') 
        if name == '':
            self.seds['default'] = opath
        else:
            self.seds[name] = opath
        return
