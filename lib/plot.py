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
from add_classes import Lightcurve, Ellipse
from functions import radio_circle, xray_circle, get_circle_size, get_symbol

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


def get_lc_dict(basepath, keys):
    lc_dict = dict()
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
    return lc_dict


def plot_fermi_flux(ax, basepath, mjd, source, **kwargs):
    lc_dict = get_lc_dict(basepath, ['ts', 'dgamma','gamma',  'flux_ul95',  'flux_err', 'flux'])
    lc_arr = dict_to_nparray(lc_dict)
    ind = np.argsort(lc_arr['tmid'])
    lc_arr = lc_arr[ind]
    mask  = (np.abs(lc_arr['ts'])>4) #& (lc_arr['flux_err'] < lc_arr['flux'])
    #gam_mask = ((lc_arr['dgamma']/lc_arr['gamma'])<2.)
    scaling = -round(np.max(np.log10(lc_arr['flux'])))+1
    ax.errorbar(lc_arr['tmid'][mask], 10**scaling*lc_arr['flux'][mask],
                 yerr = 10**scaling*lc_arr['flux_err'][mask],
                 xerr=lc_arr['bin_len'][mask]/2,linestyle=' ')
    ax.errorbar(lc_arr['tmid'][~mask], 10**scaling*lc_arr['flux_ul95'][~mask],
                 xerr=lc_arr['bin_len'][~mask]/2, color='#808080', linestyle=' ')
    serr = 0.1*(ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.errorbar(lc_arr['tmid'][~mask], 10**scaling*lc_arr['flux_ul95'][~mask],
                 yerr=serr,
                 color='#808080', uplims=True, linestyle=' ')
    ax.set_ylabel(r'$10^{'+'{:.0f}'.format(-scaling)+'}\,$'+r'ph cm$^{-2}$ s$^{-1}$')
    ax.axvline(mjd, color='#696969', linestyle='--') 
    return ax


def plot_fermi_index(ax, basepath, mjd, source, **kwargs):
    lc_dict = get_lc_dict(basepath, ['ts', 'dgamma','gamma',  'flux_ul95',  'flux_err', 'flux'])
    lc_arr = dict_to_nparray(lc_dict)
    if '4FGL' in source:
        if kwargs.get('4fgl_average', True):
            catalog = fits.open('./lib/gll_psc_v19.fit')
            names = catalog[1].data['Source_Name']
            ind = np.where(names == source)[0][0]
            av_gamma = catalog[1].data[ind]['PL_Index']
        else:
            av_gamma, _ = weighted_avg_and_std(lc_arr['gamma'],
                                               weights=1./lc_arr['dgamma'])
        ax.axhline(av_gamma,linestyle='--', color='grey',
                    alpha=0.85, zorder = -1, linewidth=0.9)
    mask  = (np.abs(lc_arr['ts'])>4) # & (lc_arr['flux_err'] < lc_arr['flux'])
    #gam_mask = ((lc_arr['dgamma']/lc_arr['gamma'])<2.)
    mask2 = (lc_arr['dgamma'] > 0)
    tot_mask = mask & mask2
    ax.errorbar(lc_arr['tmid'][mask & mask2], lc_arr['gamma'][mask & mask2],
                 yerr=lc_arr['dgamma'][mask & mask2],
                 xerr=lc_arr['bin_len'][mask& mask2]/2., linestyle='')
    ax.errorbar(lc_arr['tmid'][mask & ~mask2], lc_arr['gamma'][mask & ~mask2],
                 yerr=0, xerr=lc_arr['bin_len'][mask& ~mask2]/2., linestyle='',
                 ecolor='red')

    ax.set_ylabel('Index', labelpad=5)
    ax.set_ylim(0.0,6)
    ax.set_yticks([2,4])
    ax.axvline(mjd, color='#696969', linestyle='--')
    ax.errorbar(lc_arr['tmid'][mask & ~mask2 ], lc_arr['gamma'][mask & ~mask2],
                 yerr=0, xerr=lc_arr['bin_len'][mask& ~mask2]/2., linestyle='',
                 ecolor='red')
    return ax


def plot_radio(ax, basepath, mjd, source, **kwargs):
    idata = np.genfromtxt(basepath, delimiter=',')
    ax.errorbar(idata[:,0], idata[:,1], yerr=idata[:,2], linestyle='',
                 ecolor='red', color='red', fmt='o', ms=2)
    ax.axvline(mjd, color='#696969', linestyle='--')
    ax.set_ylabel('Jy', labelpad=5)
    return ax


def plot_xray(ax, basepath, mjd, source, **kwargs):
    idata = np.genfromtxt(basepath, delimiter=' ')
    ax.errorbar(idata[:,3], idata[:,0]*1e13, yerr=idata[:,1]*1e13, linestyle='',
                 ecolor='blue', fmt='o', color='blue',  ms=2)
    ax.axvline(mjd, color='#696969', linestyle='--')
    ax.set_ylabel(r'$10^{-13}\,$'+r'erg cm$^{-2}$ s$^{-1}$', labelpad=5)
    return ax


def make_lc_plot(lat_basepath, mjd, source, radio=None, xray=None, **kwargs):
    func_dict = {'fermi_flux': plot_fermi_flux, 'fermi_index': plot_fermi_index,
                 'radio': plot_radio, 'xray': plot_xray}

    path_dict = {'fermi_flux': lat_basepath, 'fermi_index': lat_basepath,
                 'radio': radio, 'xray': xray}
    n_panels = 2
    lcs = ['fermi_flux', 'fermi_index']
    if radio is not None:
        print('Add radio data to lightcurve')
        lcs.append('radio')
        n_panels += 1
    if xray is not None:
        print('Add x-ray data to the lightcurve')
        lcs.append('xray')
        n_panels += 1
    height_per_panel = 0.8 / n_panels
    fig = plt.figure(figsize=figsize(0.7, 0.35 * n_panels))
    ## Flux Axis
    axes = []
    xminmax = None
    for i, lc in enumerate(lcs):
        new_ax= fig.add_axes((.0, 0.2+ i*height_per_panel,1.,height_per_panel))
        new_ax= func_dict[lc](new_ax, path_dict[lc], mjd, source, **kwargs)
        axes.append(new_ax)
        if xminmax == None:
            xminmax = new_ax.get_xlim()
        else:
            new_ax.set_xlim(xminmax[0], xminmax[1])
    axes[0].set_xlabel('Time (MJD)')
    for i in range(1,len(axes)):
        axes[i].set_xticks([])
    ## Spectral Index Axis
    axes[-1].text(0.8, 1.1, source,
                  horizontalalignment='center',
                  verticalalignment='center', transform=axes[-1].transAxes)
    fig.savefig(os.path.join(lat_basepath, 'lightcurve.pdf'),
                bbox_inches='tight')
    fig.savefig(os.path.join(lat_basepath, 'lightcurve.png'),
                bbox_inches='tight', dpi=500)
    plt.close(fig)
    return


def make_sed_plot(seds_list, mw_idata=None, dec = None, twindow=None, y_min=None, y_max=None,
                  fig_scale=0.8, add_text=False, markersize=2, plot_ulims=False):
    fig, ax = newfig(fig_scale)
    ax.set_xscale('log')
    ax.set_yscale('log') 
    factor = np.log10(y_max) - np.log10(y_min)
    if mw_idata is not None:
        if len(mw_idata) > 0:
            #c = np.array(time2color(mw_idata[:,4], tmin=54500, tmax=59000))
            times = mw_idata[:,4]
            flux = mw_idata[:,1]
            flux_low = mw_idata[:,2]
            flux_up = mw_idata[:,3]
            frequency = mw_idata[:,0]
            times[times==55000] = -1
            tmask = np.array([False]*len(times))
            if twindow is not None:
                tmask = (times>twindow[0]) & (times<twindow[1])
            ulim_mask = (flux_up == flux_low)
            inds = (flux > 0) & (frequency < 1e22)
            tot_mask = inds & ~ulim_mask & ~tmask
            yerr = (flux[tot_mask] - flux_low[tot_mask],
                    flux_up[tot_mask] - flux[tot_mask])
            ax.errorbar(frequency[tot_mask] * hz_to_gev, flux[tot_mask],
                        yerr=yerr, fmt='o', color='grey', zorder=1,
                        alpha=0.5, markersize=markersize,  linestyle='')
            tot_mask = inds & ~ulim_mask & tmask
            yerr = (flux[tot_mask] - flux_low[tot_mask],
                    flux_up[tot_mask] - flux[tot_mask])
            ax.errorbar(frequency[tot_mask] * hz_to_gev, flux[tot_mask],
                        yerr=yerr, fmt='o', color='red', zorder=2,
                        alpha=0.5, markersize=markersize+0.5,  linestyle='')
            if plot_ulims:
                tot_mask = inds & ulim_mask & ~tmask
                yerr = 10**np.log10(flux[tot_mask]) - 10**(np.log10(flux[tot_mask])-0.1/8.*factor) 
                ax.errorbar(frequency[tot_mask] * hz_to_gev, flux[tot_mask],
                            yerr=yerr, fmt='o', uplims=True, color='grey', zorder=1,
                            alpha=0.5, markersize=markersize,  linestyle='')
                tot_mask = inds & ulim_mask & tmask
                yerr = 10**np.log10(flux[tot_mask]) - 10**(np.log10(flux[tot_mask])-0.1/8.*factor)   
                ax.errorbar(frequency[tot_mask] * hz_to_gev, flux[tot_mask],
                            yerr=yerr, fmt='o', uplims=True, color='red', zorder=2,
                            alpha=0.5, markersize=markersize+0.5,  linestyle='')
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
            if ts > 4 and bowtie_bool:
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
                            markersize=markersize)
                serr = 10**np.log10(yul[m]) - 10**(np.log10(yul[m])-0.1/8.*factor)
                ax.errorbar(x[m], yul[m], xerr=xerr0,
                            yerr=serr,  uplims=True,
                            color=sed_col, fmt='o', zorder=100-i,
                            markersize=markersize)


            tmin = llh['config']['selection']['tmin']
            tmax = llh['config']['selection']['tmax']
            if i == 0:
                sigma = np.max([0, pval_to_sigma(ts_to_pval(llh['sources'][source]['ts'], 1))])
                if sigma > 5:
                    sigma = np.sqrt(llh['sources'][source]['ts'])
                if add_text:
                    ax.text(0.2, 1.03, source,
                            horizontalalignment='center',
                            verticalalignment='center', transform=ax.transAxes)
                    ax.text(0.5, 1.03, '$\Sigma$: {:.1f} $\sigma$'.format(sigma),
                            horizontalalignment='center',
                            verticalalignment='center', transform=ax.transAxes)
                    ax.text(0.8, 1.03, 'MJD: {:.1f} - {:.1f}'.format(MET_to_MJD(tmin), MET_to_MJD(tmax)),
                            horizontalalignment='center',
                            verticalalignment='center', transform=ax.transAxes)
    if dec is not None:
        IC_sens = np.genfromtxt('./lib/IC_sens.txt', delimiter=',')
        inter = scipy.interpolate.UnivariateSpline(IC_sens[:,0], IC_sens[:,1], k=1,s=0)
        flux = inter(np.sin(np.radians(dec))) * MeV_to_erg * 1e6
        ax.plot([5e3, 1e6], [flux, flux], color='#0065BD', linestyle='--')
        IC_disc = np.genfromtxt('./lib/IC_disc.txt', delimiter=',')
        inter = scipy.interpolate.UnivariateSpline(IC_disc[:,0], IC_disc[:,1], k=1, s=0)
        flux = inter(np.sin(np.radians(dec))) * MeV_to_erg * 1e6
        ax.plot([5e3, 1e6], [flux, flux], color='#0065BD', linestyle='-')
    ax.set_xlim(1e-15, 1e7)
    ax.set_ylim(y_min,y_max)
    ax.set_xlabel('Energy [GeV]')
    ax.set_ylabel(r'$\nu f(\nu)$ [erg cm$^{-2}$ s$^{-1}$]')
    spath_pdf = os.path.join(seds_list[0][0], 'sed.pdf')
    spath_png = os.path.join(seds_list[0][0], 'sed.png')
    print('Save SED to {}'.format(spath_pdf))
    fig.savefig(spath_pdf, bbox_inches='tight')
    fig.savefig(spath_png, bbox_inches='tight', dpi=500)
    plt.close(fig)
    return

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


def make_ts_plot_legend(plt_basepath, srcs, max_dist):
    mask = [src.in_err for src in srcs]
    plt.clf()
    fig = plt.figure()
    fig_legend = plt.figure(figsize=(3, 2.))
    ax = fig.add_subplot(111)
    patches = []
    labels = []
    for i, src in enumerate(np.array(srcs)[mask]):
        patch = ax.plot([0], [0], linestyle='', label='{} {}'.format(i+1, src.name) ,**src.plt_style)
        labels.append('{} {}'.format(i+1, src.name))
        patches.append(patch[0])
    fig_legend.legend(patches, labels, loc='center', frameon=False, labelspacing=1.3)
    fig_legend.savefig(os.path.join(plt_basepath,'legend.png'), bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.close(fig_legend)
    return


def make_counterparts_plot(basepath, ra, dec, plt_radius, save_path='.', vou_cand=False, srcs=[], max_dist=False, legend=False, yaxis=True, error90=None):
    
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [0., 0.] ## Careful here, this is taken from the current TS definition
    wcs.wcs.cdelt = [-0.01 , 0.01]
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.ctype = ['RA---AIT', 'DEC--AIT']
    plt.clf()
    fig = plt.figure(figsize=figsize(0.5, 1.))
    ax=fig.add_axes((.0, .0,0.71,.7), projection=wcs)
    for i, src in enumerate(srcs):
        pix = get_pix_pos(wcs, src.ra, src.dec)
        if not src.in_err:
            ax.plot(pix[0], pix[1], **src.plt_style)
        else:
            ax.plot(pix[0], pix[1], label='{} {}'.format(i, src.name), **src.plt_style)
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
    ax.set_xlabel(r'R.A. (degrees)')
    if  yaxis:
        ax.set_ylabel(r'Dec. (degrees)')
    else:
        ax.set_yticklabels([])
    pix = get_pix_pos(wcs, ra, dec)
    ax.plot(pix[0] , pix[1],
            marker='o', color='k',
            ms=4, fillstyle='none', zorder=25)
    if vou_cand is not False:
        for src in vou_cand:
            pix = get_pix_pos(wcs, src[0], src[1])
            rc = radio_circle(src[2])
            xc = xray_circle(src[2])
            rcs, xcs = get_circle_size(src[2])
            if (rc is not False):
                ax.plot(pix[0], pix[1], ms=3*rcs, color=rc, marker='.')        
            
            if (xc is not False):
                ax.plot(pix[0], pix[1], ms=3*xcs,  color=xc, fillstyle='none', marker='.')
                
            symbol = get_symbol(src[2])
            if symbol is not False:
                ax.scatter(pix[0] , pix[1], **symbol)
    cpath = os.path.join(basepath, 'contour.txt')
    cpath_rad = os.path.join(basepath, 'contour_rad.txt')
    if os.path.exists(cpath):
        cdata = np.genfromtxt(cpath, delimiter=' ')
        cdata = np.vstack([cdata,cdata[0]])
        pix = get_pix_pos(wcs, np.degrees(cdata[:,0]), np.degrees(cdata[:,1]))
        ax.plot(pix[0], pix[1], color='b', linewidth=0.5)
    elif os.path.exists(cpath_rad):
        print('Use dfit Contour')
        cdata = np.genfromtxt(cpath_rad, delimiter=' ', skip_header=1)
        cdata = np.vstack([cdata,cdata[0]])
        print np.degrees(cdata[:,0])
        print 360+np.degrees(cdata[:,1])
        pix = get_pix_pos(wcs,  360+np.degrees(cdata[:,1]), np.degrees(cdata[:,0]))
        ax.plot(pix[0], pix[1], color='b', linewidth=0.5)
    else:
        if error90 is not None:
            if isinstance(error90, float):
                vertex = (ra*u.degree, dec*u.degree) #long, lat
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
        
    if legend:
        ax.legend(bbox_to_anchor=(1.1, 0.5), loc='center left', prop={'size': 9})
    ax.grid(color='k', ls='--')
    factorm = 1./wcs.wcs.cdelt[0]
    factorp = 1./wcs.wcs.cdelt[1]
    print plt_radius*factorm, plt_radius*factorp
    ax.set_xlim(plt_radius*factorm, plt_radius*factorp)
    ax.set_ylim(plt_radius*factorm, plt_radius*factorp)    
    print('Save Counterpart Plot')
    fig.savefig(os.path.join(save_path, 'counterparts.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(save_path, 'counterparts.png'), bbox_inches='tight', dpi=300)
    return




def make_ts_plot(plt_basepath, srcs, vou_cand, bf_ra=None, bf_dec=None, plt_mode='tsmap', legend=False, yaxis=True, error90=None):
    plt.clf()
    fname = 'fit1_pointsource_powerlaw_2.00_{}.fits'.format(plt_mode)
    fits_path = os.path.join(plt_basepath, fname)
    if os.path.exists(fits_path):
        inp = fits.open(fits_path)
    else:
        print('{} not yet ready'.format(fits_path))
        return
    wcs = WCS(inp[2].header)
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
    if bf_ra is None:
        bf_ra = hdu.header['CRVAL1']
        bf_dec = hdu.header['CRVAL2']
    fig = plt.figure(figsize=figsize(0.4, 1.))
    ax=fig.add_axes((.0, .0,0.7,.7), projection=wcs)
    cpath = os.path.join(plt_basepath, '../contour.txt')
    cpath_rad = os.path.join(plt_basepath, '../contour_df.txt')
    if os.path.exists(cpath):
        cdata = np.genfromtxt(cpath, delimiter=' ')
        cdata = np.vstack([cdata,cdata[0]])
        pix = get_pix_pos(wcs, np.degrees(cdata[:,0]), np.degrees(cdata[:,1]))
        ax.plot(pix[0], pix[1], color='b', linewidth=0.5)
    elif os.path.exists(cpath_rad):
        print('Use dfit Contour')
        cdata = np.genfromtxt(cpath_rad, delimiter=' ', skip_header=1)
        cdata = np.vstack([cdata,cdata[0]])
        pix = get_pix_pos(wcs, 360 + np.degrees(cdata[:,1]), np.degrees(cdata[:,0]))
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
    if cpath: #os.path.exists(cpath):
        cdata = np.genfromtxt(cpath, delimiter=' ')
        cdata = np.vstack([cdata,cdata[0]])
        pix = get_pix_pos(wcs, cdata[:,0], cdata[:,1])
        ax.plot(pix[0], pix[1], color='k', linewidth=0.5)
    for i, src in enumerate(srcs):
        pix = get_pix_pos(wcs, src.ra, src.dec)
        ax.plot(pix[0], pix[1],label=src.name, **src.plt_style)
    if vou_cand is not False:
        for src in vou_cand:
            pix = get_pix_pos(wcs, src[0], src[1])
            rc = radio_circle(src[2])
            xc = xray_circle(src[2])
            rcs, xcs = get_circle_size(src[2])
            if (rc is not False):
                ax.plot(pix[0], pix[1], ms=3*rcs, color=rc, marker='.')

            if (xc is not False):
                ax.plot(pix[0], pix[1], ms=3*xcs,  color=xc, fillstyle='none', marker='.')

            symbol = get_symbol(src[2])
            if symbol is not False:
                ax.plot(pix[0] , pix[1], **symbol)
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
    levels=np.array([2,3,4])
    CS1 = ax.contour(Z, levels=levels,
                    colors='black', linewidths=(0.2,))
    levels=np.array([5])
    if np.max(Z) > 5:
        CS2 = ax.contour(Z, levels=levels,
                        colors='black', linewidths=(0.8,))
    ax.set_xlabel(r'R.A. (degrees)')
    if  yaxis:
        ax.set_ylabel(r'Dec. (degrees)')
    else:
        ax.set_yticklabels([])
    pix = get_pix_pos(wcs, bf_ra, bf_dec)
    ax.plot(pix[0] , pix[1],
            marker='o', color='blue',
            ms=3, fillstyle='none', zorder=2)
    
    ax.set_xlim(0, hdu.data.shape[1])
    ax.set_ylim(0, hdu.data.shape[0]) 
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
  #  plt.tight_layout()
    fig.savefig(os.path.join(plt_basepath,'{}.png'.format(plt_mode)),
                bbox_inches='tight', dpi=500)
    plt.close(fig)
    return


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
