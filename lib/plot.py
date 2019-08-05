# coding: utf-8

import sys
sys.path.append('/home/ga53lag/Software/python_scripts/')
from fancy_plot import *
import numpy as np
import os
from myfunctions import dict_to_nparray, MET_to_MJD, GreatCircleDistance, pval_to_sigma, ts_to_pval
from astropy.time import Time
import pyfits as fits
import warnings
from matplotlib.colors import LinearSegmentedColormap
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.modeling.functional_models import Ellipse2D
import astropy.units as u
from collections import OrderedDict
import scipy.interpolate
from astropy.wcs import WCS
import yaml
from regions import EllipseSkyRegion
from astropy.coordinates import SkyCoord
from source_class import Ellipse


markers = ['o', 's', 'P', 'p', '*' , 'x', 'X', 'D', 4, 5, 6, 7, 'H','d', 'v' ,'^', '<', '>', 1, 2, 3 ,8, '+' ,'h']
ts_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', '#800000'])
re_cmap = LinearSegmentedColormap.from_list('mycmap2', ['#67a9cf', '#f7f7f7', '#ef8a62'])
mw_map= LinearSegmentedColormap.from_list('mycmap3', ['#bdbdbd','#939393' ,'red'])
hz_to_gev = 4.135*10**(-24)
MeV_to_erg = 1.60218e-6
ul_ts_threshold = 4


def time2color(ts, tmin = -1, tmax = -1):
    if tmin == -1:
        tmin = np.min(ts[ts>0])
    if tmax == -1:
        tmax = np.max(ts)
    delta_t = tmax - tmin
    c = []
    for t in ts:
        if t<tmin:
            col = 0.
        elif t> tmax:
            col =1.
        else:
            col = ((t-tmin)/delta_t)
        c.append(mw_map(col))
    return c


def make_lc_plot(basepath, mjd, **kwargs):
    lc_dict = dict()
    keys = ['ts', 'dgamma','gamma',  'flux_ul95',  'flux_err', 'flux']
    for key in keys:
        lc_dict[key] = [] 
    lc_dict['tmid'] = []
    lc_dict['bin_len'] = []
    for folder in os.listdir(basepath):
        path = os.path.join(basepath,folder,'llh.npy')
        if os.path.exists(path):
            inp = np.load(path, allow_pickle=True)[()]
            source = inp['config']['selection']['target']
            flux_dict = inp['sources'][source]
        else:
            continue
        if not np.isfinite(flux_dict['ts']):
            continue
        for key in keys:
            if key=='gamma':
                lc_dict[key].append(get_index(flux_dict))
            elif key=='dgamma':
                lc_dict[key].append(get_index_err(flux_dict))   
            else:
                lc_dict[key].append(flux_dict[key])
        lc_dict['tmid'].append((float(folder.split('_')[1])+float(folder.split('_')[0]))/2)
        lc_dict['bin_len'].append(float(folder.split('_')[1])-float(folder.split('_')[0]))
    lc_arr = dict_to_nparray(lc_dict)
    ind = np.argsort(lc_arr['tmid'])
    lc_arr = lc_arr[ind]
    fig = plt.figure(figsize=figsize(0.5, 0.7))
    mask  = (np.abs(lc_arr['ts'])>4) & (lc_arr['flux_err'] < lc_arr['flux'])
    gam_mask = ((lc_arr['dgamma']/lc_arr['gamma'])<2.)
    ## Flux Axis
    ax1=fig.add_axes((.0, .20,1.,.4))
    scaling = -round(np.max(np.log10(lc_arr['flux'])))+1
    ax1.errorbar(lc_arr['tmid'][mask], 10**scaling*lc_arr['flux'][mask],
                 yerr = 10**scaling*lc_arr['flux_err'][mask],
                 xerr=lc_arr['bin_len'][mask]/2,linestyle=' ')
    ax1.errorbar(lc_arr['tmid'][~mask], 10**scaling*lc_arr['flux_ul95'][~mask],
                 xerr=lc_arr['bin_len'][~mask]/2, color='#808080', linestyle=' ')
    serr = 0.1*(ax1.get_ylim()[1] - ax1.get_ylim()[0])
    ax1.errorbar(lc_arr['tmid'][~mask], 10**scaling*lc_arr['flux_ul95'][~mask],
                 yerr=serr,
                 color='#808080', uplims=True, linestyle=' ')
    ax1.set_ylabel(r'$10^{'+'{:.0f}'.format(-scaling)+'}\,$'+r'ph cm$^{-2}$ s$^{-1}$')
    ax1.axvline(mjd, color='#696969', linestyle='--')
    ax1.set_xlabel('Time (MJD)')

    ## Spectral Index Axis
    ax2=fig.add_axes((.0, .6,1.,.4))
    if '4FGL' in source:
        if kwargs.get('4fgl_average', True):
            catalog = fits.open('./lib/gll_psc_v19.fit')
            names = catalog[1].data['Source_Name']
            ind = np.where(names == source)[0][0]
            av_gamma = catalog[1].data[ind]['PL_Index']
        else:
            av_gamma, _ = weighted_avg_and_std(lc_arr['gamma'],
                                               weights=1./lc_arr['dgamma'])
        ax2.axhline(av_gamma,linestyle='--', color='grey',
                    alpha=0.85, zorder = -1, linewidth=0.9)

    mask2 = (lc_arr['dgamma'] > 0)
    ax2.errorbar(lc_arr['tmid'][mask & mask2 & gam_mask], lc_arr['gamma'][mask & mask2 & gam_mask],
                 yerr=lc_arr['dgamma'][mask & mask2 & gam_mask],
                 xerr=lc_arr['bin_len'][mask& mask2 &gam_mask]/2., linestyle='')
    ax2.errorbar(lc_arr['tmid'][mask & ~mask2 & gam_mask], lc_arr['gamma'][mask & ~mask2 & gam_mask],
                 yerr=0, xerr=lc_arr['bin_len'][mask& ~mask2 &gam_mask]/2., linestyle='',
                 ecolor='red')
    ax2.set_ylabel('Index', labelpad=5)
    ax2.set_ylim(0.0,6)
    ax2.set_xticks([])
    ax2.set_yticks([2,4])
    ax2.axvline(mjd, color='#696969', linestyle='--')
    ax2.set_xlim(ax1.get_xlim()[0], ax1.get_xlim()[1])
    ax2.text(0.8, 1.1, source,
        horizontalalignment='center',
        verticalalignment='center', transform=ax2.transAxes)
    plt.savefig(os.path.join(basepath, 'lightcurve.pdf'),
                bbox_inches='tight')
    return


def make_sed_plot(seds_list, mw_data=None, dec = None, twindow=None):
    fig, ax = newfig(0.9)
    ax.set_xscale('log')
    ax.set_yscale('log') 
    y_vals = []
    if mw_data is not None:
        try:
            mw_idata = np.atleast_2d(np.genfromtxt(mw_data, skip_header=1, usecols=(0,1,2,3,4)))
            if len(mw_idata) > 0:
                inds = (mw_idata[:,1] > 0) & (mw_idata[:,0] < 1e22)
                y_vals.extend(mw_idata[:,1][inds])
        except Exception as inst:
            pass
    for i, sed_list in enumerate(seds_list):
        basepath = sed_list[0]
        if os.path.exists(os.path.join(basepath, 'sed.npy')):
            sed = np.load(os.path.join(basepath, 'sed.npy'), allow_pickle=True)[()]
        else:
            continue
        m = sed['ts'] < ul_ts_threshold
        y_vals.extend(sed['e2dnde'][~m])
        y_vals.extend(sed['e2dnde_ul95'][m])
    if len(y_vals) ==0:
        factor = 1.
    else:
        factor = np.log10(np.max(y_vals)) - np.log10(np.min(y_vals))
    if mw_data is not None:
        try:
            mw_idata = np.atleast_2d(np.genfromtxt(mw_data, skip_header=1, usecols=(0,1,2,3,4)))
        except Exception as inst:
            mw_idata = np.array([])
        if len(mw_idata) > 0:
            #c = np.array(time2color(mw_idata[:,4], tmin=54500, tmax=59000))
            times = mw_idata[:,4]
            tmask = np.array([False]*len(times))
            if twindow is not None:
                tmask = (times>twindow[0]) & (times<twindow[1])
            ulim_mask = (mw_idata[:,2] == mw_idata[:,3])
            inds = (mw_idata[:,1] > 0) & (mw_idata[:,0] < 1e22)
            tot_mask = inds & ~ulim_mask & ~tmask
            yerr = (mw_idata[:,1][tot_mask] - mw_idata[:,2][tot_mask],
                    mw_idata[:,3][tot_mask] - mw_idata[:,1][tot_mask])
            ax.errorbar(mw_idata[:,0][tot_mask] * hz_to_gev, mw_idata[:,1][tot_mask],
                        yerr=yerr, fmt='o', color='grey', zorder=1,
                        alpha=0.5, markersize=3,  linestyle='')
            tot_mask = inds & ~ulim_mask & tmask
            yerr = (mw_idata[:,1][tot_mask] - mw_idata[:,2][tot_mask],
                    mw_idata[:,3][tot_mask] - mw_idata[:,1][tot_mask])
            ax.errorbar(mw_idata[:,0][tot_mask] * hz_to_gev, mw_idata[:,1][tot_mask],
                        yerr=yerr, fmt='o', color='red', zorder=2,
                        alpha=0.5, markersize=3,  linestyle='')

            tot_mask = inds & ulim_mask & ~tmask
            yerr = 10**np.log10(mw_idata[:,1][tot_mask]) - 10**(np.log10(mw_idata[:,1][tot_mask])-0.1/8.*factor) 
            ax.errorbar(mw_idata[:,0][tot_mask] * hz_to_gev, mw_idata[:,1][tot_mask],
                        yerr=yerr, fmt='o', uplims=True, color='grey', zorder=1,
                        alpha=0.5, markersize=3,  linestyle='')
            tot_mask = inds & ulim_mask & tmask
            yerr = 10**np.log10(mw_idata[:,1][tot_mask]) - 10**(np.log10(mw_idata[:,1][tot_mask])-0.1/8.*factor)   
            ax.errorbar(mw_idata[:,0][tot_mask] * hz_to_gev, mw_idata[:,1][tot_mask],
                        yerr=yerr, fmt='o', uplims=True, color='red', zorder=2,
                        alpha=0.5, markersize=3,  linestyle='')
     
    for i, sed_list in enumerate(seds_list):
        basepath = sed_list[0]
        sed_col = sed_list[1]
        bowtie_col = sed_list[2]
        bowtie_fill = sed_list[3]
        bowtie_bool = sed_list[4]
        sed_bool = sed_list[5]
        if os.path.exists(os.path.join(basepath, 'sed.npy')) & \
           os.path.exists(os.path.join(basepath, 'bowtie.npy')) & \
           os.path.exists(os.path.join(basepath, 'llh.npy')): 
            sed = np.load(os.path.join(basepath, 'sed.npy'), allow_pickle=True)[()]
            bowtie = np.load(os.path.join(basepath, 'bowtie.npy'),allow_pickle=True)[()]
            llh = np.load(os.path.join(basepath, 'llh.npy'), allow_pickle=True)[()]
            source = llh['config']['selection']['target']
            ts = llh['sources'][source]['ts']
            e2 = 10 ** (2 * bowtie['log_energies'])
            energies = np.array(10 ** bowtie['log_energies'])/1e3
            m = sed['ts'] < ul_ts_threshold
            x = sed['e_ctr']/1e3
            y = sed['e2dnde']*MeV_to_erg
            yerr = sed['e2dnde_err']*MeV_to_erg
            yerr_lo = sed['e2dnde_err_lo']*MeV_to_erg
            yerr_hi = sed['e2dnde_err_hi']*MeV_to_erg
            yul = sed['e2dnde_ul95']*MeV_to_erg
            delo = (sed['e_ctr'] - sed['e_min'])/1e3
            dehi = (sed['e_max'] - sed['e_ctr'])/1e3
            xerr0 = np.vstack((delo[m], dehi[m]))
            xerr1 = np.vstack((delo[~m], dehi[~m]))
            if ts > 9 and bowtie_bool:
                ax.plot(energies,
                        bowtie['dnde'] * e2 * MeV_to_erg, color=bowtie_col,
                        zorder=10-i)
                ax.plot(energies,
                        bowtie['dnde_lo'] * e2 * MeV_to_erg, color=bowtie_col,
                        linestyle='--', zorder=10-i)
                ax.plot(energies,
                        bowtie['dnde_hi'] * e2 * MeV_to_erg, color=bowtie_col,
                        linestyle='--', zorder=10-i)
                if bowtie_fill:
                    ax.fill_between(energies,
                                    bowtie['dnde_lo'] * e2 * MeV_to_erg,
                                    bowtie['dnde_hi'] * e2 * MeV_to_erg,
                                    alpha=0.5, color=bowtie_col, zorder=10-i)
            if sed_bool:
                ax.errorbar(x[~m], y[~m], xerr=xerr1,
                            yerr=(yerr_lo[~m], yerr_hi[~m]),
                            linestyle='', color=sed_col, fmt='o', zorder=100-i,
                            markersize=3)
                serr = 10**np.log10(yul[m]) - 10**(np.log10(yul[m])-0.1/8.*factor)
                ax.errorbar(x[m], yul[m], xerr=xerr0,
                            yerr=serr,  uplims=True,
                            color=sed_col, fmt='o', zorder=100-i,
                            markersize=3)


            tmin = llh['config']['selection']['tmin']
            tmax = llh['config']['selection']['tmax']
            if i == 0:
                sigma = np.max([0, pval_to_sigma(ts_to_pval(llh['sources'][source]['ts'], 1))])
                if sigma > 5:
                    sigma = np.sqrt(llh['sources'][source]['ts'])
                ax.text(0.2, 1.02, source,
                        horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
                ax.text(0.5, 1.02, '$\Sigma$: {:.1f} $\sigma$'.format(sigma),
                        horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
                ax.text(0.8, 1.02, 'MJD: {:.1f} - {:.1f}'.format(MET_to_MJD(tmin), MET_to_MJD(tmax)),
                        horizontalalignment='center',
                        verticalalignment='center', transform=ax.transAxes)
    if dec is not None:
        IC_sens = np.genfromtxt('./lib/IC_sens.txt', delimiter=',')
        inter = scipy.interpolate.UnivariateSpline(IC_sens[:,0], IC_sens[:,1], k=1,s=0)
        flux = inter(np.sin(np.radians(dec))) * MeV_to_erg * 1e6
        ax.plot([5e3, 1e6], [flux, flux], color='green', linestyle='--')
        IC_disc = np.genfromtxt('./lib/IC_disc.txt', delimiter=',')
        inter = scipy.interpolate.UnivariateSpline(IC_disc[:,0], IC_disc[:,1], k=1, s=0)
        flux = inter(np.sin(np.radians(dec))) * MeV_to_erg * 1e6
        ax.plot([5e3, 1e6], [flux, flux], color='green', linestyle='-')
    if np.min(y_vals) < 1e-17:
        print('Warning: Minimum data point could be out of the plotting range ({})'.format(np.min(y_vals)))
    ax.set_xlim(1e-15, 1e7)
    ax.set_ylim(np.min(y_vals),1e-9)
    ax.set_xlabel('Energy [GeV]')
    ax.set_ylabel(r'$\nu f(\nu)$ [erg cm$^{-2}$ s$^{-1}$]')
    plt.tight_layout()
    spath = os.path.join(seds_list[0][0], 'sed.pdf')
    print('Save SED to {}'.format(spath))
    fig.savefig(spath, bbox_inches='tight')


def make_edges(data_f):
    header = data_f[0].header
    xmin = header['CRVAL1']-header['CRPIX1']*abs(header['CDELT1'])
    xmax = header['CRVAL1']+header['CRPIX1']*abs(header['CDELT1'])
    ymin = header['CRVAL2']-header['CRPIX2']*abs(header['CDELT2'])
    ymax = header['CRVAL2']+header['CRPIX2']*abs(header['CDELT2']) 
    xbins = header['NAXIS1']
    ybins = header['NAXIS2']
    return xmin, xmax, ymin, ymax, xbins, ybins


def get_pix_pos(wcs, ra, dec):
    pix = wcs.wcs_world2pix(np.array(zip(np.atleast_1d(ra), np.atleast_1d(dec))),1)
    return pix[:,0]-1, pix[:,1]-1


def make_ts_plot_legend(plt_basepath, srcs):
    fig = plt.figure()
    fig_legend = plt.figure(figsize=(2, 1.25))
    ax = fig.add_subplot(111)
    patches = []
    labels = []
    for i, src in enumerate(srcs):
        patch = ax.scatter([0], [0],  marker=markers[i],
                           label=src.name, color='k')
        labels.append(src.name)
        patches.append(patch)
    fig_legend.legend(patches, labels, loc='center', frameon=False)
    fig_legend.savefig(os.path.join(plt_basepath,'legend.png'), bbox_inches='tight', dpi=300)
    return

def make_ts_plot(plt_basepath, srcs, vou_cand, plt_mode='tsmap', legend=False, yaxis=True, error90=None):
    fname = 'fit1_pointsource_powerlaw_2.00_{}.fits'.format(plt_mode)
    fits_path = os.path.join(plt_basepath, fname)
    if os.path.exists(fits_path):
        inp = fits.open(fits_path)
        try:
            with open(os.path.join(plt_basepath, 'config.yaml'), 'r') as f:
                yml = yaml.safe_load(f)
            binning = yml['binning']
            delta_pix = (binning[ 'roiwidth']/binning['binsz'] - 6./binning['binsz']) / 2. ## Assume radius of 3 deg (total width 6 deg)
        except Exception as inst:
            delta_pix = 10 # random default value
            print inst
    else:
        print('{} not yet ready'.format(fits_path))
        return
    wcs = WCS(inp[2].header)
    cand_pos = np.genfromtxt(vou_cand)
    cand = {'ra': cand_pos[:,0], 'dec': cand_pos[:,1]}
    if (len(cand['ra'])>0) and (len(srcs) >0):
        distances = np.ones(len(cand['ra']))
        for i in range(len(cand['ra'])):
            tdist= [GreatCircleDistance(
                      cand['ra'][i], cand['dec'][i], srcs[j].ra,
                      srcs[j].dec, unit='deg') for j in range(len(srcs))]
            distances[i] = np.degrees(np.min(tdist))
        mask = distances>0.1
        inds =  np.argsort(distances)[len(srcs):]
        cand['ra'] = cand['ra'][mask]
        cand['dec'] =cand['dec'][mask]
 
    if plt_mode == 'tsmap':
        hdu = inp[2]
        Z=pval_to_sigma(ts_to_pval(hdu.data,1.))
        minmax = (0,8)
        ticks = 2
        cmap = ts_cmap
    elif plt_mode == 'residmap':
        hdu = inp[0]
        Z=hdu.data
        minmax = (-6, 6)
        ticks = 2 
        cmap = re_cmap

    bfpath = os.path.join(plt_basepath, '../bf.txt')
    if os.path.exists(bfpath):
        bfdata = np.genfromtxt(bfpath, delimiter=',')
        bf_ra = bfdata[0]
        bf_dec = bfdata[1]
    else:
        bf_ra = hdu.header['CRVAL1']
        bf_dec = hdu.header['CRVAL2']
    fig = plt.figure(figsize=figsize(0.4, 1.))
    plt.clf()
    ax=fig.add_axes((.0, .0,0.7,.7), projection=wcs)
    cpath = os.path.join(plt_basepath, '../contour.txt')
    if os.path.exists(cpath):
        cdata = np.genfromtxt(cpath, delimiter=',')
        cdata = np.vstack([cdata,cdata[0]])
        pix = get_pix_pos(wcs, cdata[:,0], cdata[:,1])
        ax.plot(pix[0], pix[1], color='b', linewidth=0.5)
    elif error90 is not None:
        if isinstance(error90, float):
            vertex = (hdu.header['CRVAL1']*u.degree,hdu.header['CRVAL2']*u.degree) #long, lat
            x =  SphericalCircle(vertex, error90*u.degree) 
            cdata = x.get_xy()
            pix = get_pix_pos(wcs, cdata[:,0], cdata[:,1])
            ax.plot(pix[0], pix[1], color='b', linewidth=0.5)
        elif isinstance(error90, Ellipse):
            ell= EllipseSkyRegion(SkyCoord(error90.center_ra * u.deg, error90.center_dec * u.deg, frame='icrs'),
                                  2. * error90.ra_ax * u.deg, 2. * error90.dec_ax * u.deg,
                                  angle = (error90.rotation) * u.deg )
            pix = ell.to_pixel(wcs)
            pix.plot(ax=ax, color='b')
    else:
        vertex = (hdu.header['CRVAL1']*u.degree,hdu.header['CRVAL2']*u.degree) #long, lat
        x =  SphericalCircle(vertex,1.5*u.degree) 
        cdata = x.get_xy()
        pix = get_pix_pos(wcs, cdata[:,0], cdata[:,1])
        ax.plot(pix[0], pix[1], color='b', linewidth=0.5)
        print('Use new circle')

    cpath = os.path.join(plt_basepath, '../contour50.txt')
    if os.path.exists(cpath):
        cdata = np.genfromtxt(cpath, delimiter=',')
        cdata = np.vstack([cdata,cdata[0]])
        pix = get_pix_pos(wcs, cdata[:,0], cdata[:,1])
        ax.plot(pix[0], pix[1], color='k', linewidth=0.5)
    for i, src in enumerate(srcs):
        pix = get_pix_pos(wcs, src.ra, src.dec)
        ax.plot(pix[0], pix[1],
                marker=markers[i],
                linestyle = '',
                label=src.name,
                color='k', ms=4, zorder=3)
    if len(cand['ra'])>0:
        pix = get_pix_pos(wcs, cand['ra'], cand['dec'])
        ax.plot(pix[0], pix[1], color='#a8a8a8', ms=4, marker="8", zorder=1,
                linestyle = '', fillstyle='none', label='VOU Sources', mew=1)
    lon = ax.coords[0]
    lat = ax.coords[1]
    lat.set_ticks_position('l')
    lat.set_ticklabel_position('l')
    lat.set_axislabel_position('l')
    lon.set_ticks_position('bt')
    lon.set_ticklabel_position('bt')
    lon.set_axislabel_position('b')
    lon.set_major_formatter('d.d')
    lat.set_major_formatter('d.d') 
    cbar = ax.contourf(Z, levels=np.linspace(minmax[0], minmax[1], 500),
                       cmap=cmap)
    levels=np.array([2,3,4,5])
    CS = ax.contour(Z, levels=levels,
                    colors='black', linewidths=(0.3,))
    ax.set_xlabel(r'R.A. (degrees)')
    if  yaxis:
        ax.set_ylabel(r'Dec. (degrees)')
    else:
        ax.set_yticklabels([])
    pix = get_pix_pos(wcs, bf_ra, bf_dec)
    ax.plot(pix[0] , pix[1],
            marker='o', color='blue',
            ms=3, fillstyle='none', zorder=2)
    
    #ax.set_xlim(0, hdu.data.shape[1])
    #ax.set_ylim(0, hdu.data.shape[0]) 
    ax2=fig.add_axes((.0, .80,0.7, 0.05))
    plt_cbar = fig.colorbar(cbar, orientation="horizontal", cax=ax2,
                            ticks=np.arange(minmax[0], minmax[1] , ticks))
    if plt_mode == 'tsmap':
        plt_cbar.set_label(r'Significance [$\sigma$]', labelpad=8)
    elif plt_mode == 'residmap':
        plt_cbar.set_label(r'Significance [$\sigma$]', labelpad=8)
    plt_cbar.ax.xaxis.set_ticks_position('top')
    plt_cbar.ax.xaxis.set_label_position('top')
    ax2.tick_params(axis='x', which='major', pad=3)
    if legend:
        ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', prop={'size': 9})

    ax.grid(color='k', ls='--')
    plt.tight_layout()
    plt.savefig(os.path.join(plt_basepath,'{}.png'.format(plt_mode)),
                bbox_inches='tight', dpi=300)


def get_index(flux_dict):
    if flux_dict['SpectrumType']=='PowerLaw':
        return flux_dict['spectral_pars']['Index']['value']
    if flux_dict['SpectrumType']=='LogParabola':
        return flux_dict['spectral_pars']['alpha']['value'] 

def get_index_err(flux_dict):
    if flux_dict['SpectrumType']=='PowerLaw':
        ret = flux_dict['spectral_pars']['Index']['error']
    elif flux_dict['SpectrumType']=='LogParabola':
        ret = flux_dict['spectral_pars']['alpha']['error']
    if not np.isfinite(ret):
        ret = -1
    return ret

 
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))
