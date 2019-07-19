import numpy as np
from roi_functions import GreatCircleDistance, vou_path
from astropy.io import fits

class Source(object):
    """Class holding all the information for a Counterpart Candidate"""

    def __init__(self, name, ra, dec):
        self.name = name
        self.names = []
        self.ra = ra  
        self.dec = dec
        self.lightcurves = []
        
    def make_sed_plot(self):
        # To be implemented
        return 

    def make_lc_plot(self):
        # To be implemented
        return 
    
    def make_fixed_binning_lightcurve(self, threshold=1e3):
        # To be implemented
        return
    
    def make_sed(self):
        # To be implemented
        return
    
    def make_summary(self):
        # To be implemented
        return 
    
    def pick_primary_name(self):
        # To be implemented
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
        if not os.path.exists(self.path):
            os.makedirs(self.bpath)
        self.sed_path = os.path.join(bpath_src, path_settings['sed'])
        if not os.path.exists(self.sed_path):
            os.makedirs(self.sed_path)
        self.lc_path = os.path.join(bpath_src, path_settings['lc'])
        if not os.path.exists(self.lc_path):
            os.makedirs(self.lc_path)
        self.mw_data_path = os.path.join(bpath_src, 'sed.txt')
        if not os.path.exists(self.mw_data_path):
            os.makedirs(self.mw_data_path)
        return

    def get_mw_data(self):
        if '4FGL' in src.name:
            data = fits.open('gll_psc_v19.fit')[1].data
            ind = np.where(data['Source_Name']==src.name)[0]
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
                os.system('{vou_path} {ra} {dec} {loc_str}'.format(vou_path=vou_path, ra=src.ra,
                                                                   dec=src.dec,loc_str=loc_str))
            if not os.path.exists('candidates.eps'):
                fname_new = os.path.join(ofolder, src.name.replace(' ', '_').replace('.','_') + '.ps')
                cand_ps = os.path.join(ofolder, 'candidates.ps')
                os.rename(cand_ps, fname_new)
                os.system('ps2eps -B ' + fname_new)
            else:
                cand_eps = os.path.join(ofolder, 'candidates.eps')
                fname_new = os.path.join(ofolder, src.name.replace(' ', '_').replace('.','_') + '.eps')
                os.rename(cand_eps, fname_new)
            os.chdir(this_path)
            cp_info = np.array(np.genfromtxt(os.path.join(ofolder, 'candidates_posix.txt'), dtype=float), ndmin=2)
            if len(cp_info) == 0:
                cp_ra = src.ra
                cp_dec = src.dec
            else:
                try:
                    cp_ra = cp_info[:,0][0]
                    cp_dec = cp_info[:,1][0]
                    print('Found Counterpart at ra {}, dec {}'.format(cp_ra, cp_dec))
                except Exception as inst:
                    cp_ra = src.ra
                    cp_dec = src.dec
        else:
            cp_ra = src.ra
            cp_dec = src.dec
        pre_str = '{vou_path} {ra} {dec} {loc_str} -s ; cat Sed.txt > {bpath}'
        for i in range(2):
            os.system(pre_str.format(vou_path=vou_path, ra=cp_ra, dec=cp_dec,  bpath=os.path.join(bpath_src,'sed.txt'), loc_str=2))
        return
