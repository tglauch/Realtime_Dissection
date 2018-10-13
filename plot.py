# coding: utf-8

import sys
sys.path.append('/home/ga53lag/Software/python_scripts/')
from fancy_plot import *
import numpy as np
import os
from myfunctions import dict_to_nparray, MET_to_MJD
from astropy.time import Time
import pyfits as fits
from matplotlib.colors import LinearSegmentedColormap
ts_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', '#800000'])
re_cmap = LinearSegmentedColormap.from_list('mycmap2', ['#67a9cf', '#f7f7f7', '#ef8a62'])
MeV_to_erg = 1.60218e-6

def make_lc_plot(basepath):
    lc_dict = dict()
    keys = ['dgamma','gamma', 'ts',  'flux_ul95',  'flux_err', 'flux']
    for key in keys:
        lc_dict[key] = [] 
    lc_dict['tmid'] = []
    lc_dict['bin_len'] = []
    for folder in os.listdir(basepath):
        path = os.path.join(basepath,folder,'llh.npy')
        if os.path.exists(path):
            inp = np.load(path)[()]
            source = inp['config']['selection']['target']
            print(source)
            flux_dict = inp['sources'][source]
        else:
            continue
        lc_dict['tmid'].append((float(folder.split('_')[1])+float(folder.split('_')[0]))/2)
        lc_dict['bin_len'].append(float(folder.split('_')[1])-float(folder.split('_')[0]))
        for key in keys:
            if key=='gamma':
                lc_dict[key].append(get_index(flux_dict))
            elif key=='dgamma':
                lc_dict[key].append(get_index_err(flux_dict))   
            else:
                lc_dict[key].append(flux_dict[key])
    lc_arr = dict_to_nparray(lc_dict)
    ind = np.argsort(lc_arr['tmid'])
    lc_arr = lc_arr[ind]
    fig = plt.figure(figsize=figsize(0.5, 0.7))
    mask  = (np.abs(lc_dict['ts'])>4)

    ## Flux Axis
    ax1=fig.add_axes((.0, .60,1.,.4))
    scaling = -round(np.max(np.log10(lc_arr['flux'])))+1
    ax1.errorbar(lc_arr['tmid'], 10**scaling*lc_arr['flux'],
                 yerr = 10**scaling*lc_arr['flux_err'],
                 xerr=lc_arr['bin_len']/2,linestyle=' ')
    ax1.errorbar(lc_arr['tmid'][~mask], lc_arr['flux_ul95'][~mask],
                 xerr=lc_arr['bin_len'][~mask]/2, color='#808080')
    ax1.errorbar(lc_arr['tmid'][~mask], lc_arr['flux_ul95'][~mask],
                 yerr=0.1*np.ones(len(lc_arr['flux_ul95']))[~mask],
                 color='#808080', uplims=True)
    ax1.set_ylabel(r'$10^{'+'{:.0f}'.format(-scaling)+'}\,$'+r'ph cm$^{-2}$ s$^{-1}$')
    ax1.set_xticks([])

    ## Spectral Index Axis
    ax2=fig.add_axes((.0, .2,1.,.40))
    av_gamma, std_gamma = weighted_avg_and_std(lc_arr['gamma'],
                                               weights=1./lc_arr['dgamma'])
    ax2.axhline(av_gamma,linestyle='--', color='grey',
                alpha=0.85, zorder = -1, linewidth=0.9)
    ax2.errorbar(lc_arr['tmid'][mask], lc_arr['gamma'][mask],
                 yerr=lc_arr['dgamma'][mask],
                 xerr=lc_arr['bin_len'][mask]/2., linestyle='')
    ax2.set_ylabel('Index', labelpad=5)
    ax2.set_xlabel('Time (MJD)')
    ax1.text(0.8, 1.1, source,
        horizontalalignment='center',
        verticalalignment='center', transform=ax1.transAxes)
    plt.savefig(os.path.join(basepath, 'lightcurve.pdf'),
                bbox_inches='tight')
    return

def make_sed_plot(basepath):
    sed = np.load(os.path.join(basepath, 'sed.npy'))[()]
    bowtie = np.load(os.path.join(basepath, 'bowtie.npy'))[()]
    llh = np.load(os.path.join(basepath, 'llh.npy'))[()]
    fig, ax = newfig(0.9)
    ax = fig.add_subplot(111)
    color='red'
    e2 = 10 ** (2 * bowtie['log_energies'])
    energies = np.array(10 ** bowtie['log_energies'])/1e3
    ul_ts_threshold = 4
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

    ax.errorbar(x[~m], y[~m], xerr=xerr1,
                 yerr=(yerr_lo[~m], yerr_hi[~m]),
                 linestyle='', color='k', fmt='o', zorder=2)
    ax.errorbar(x[m], yul[m], xerr=xerr0,
                 yerr=(yul[m]*(0.15), np.zeros(len(yul[m]))),  uplims=True,
                 color='k', fmt='o', zorder=2)
    ax.plot(energies,
            bowtie['dnde'] * e2 * MeV_to_erg, color=color, zorder=1)

    ax.plot(energies,
            bowtie['dnde_lo'] * e2 * MeV_to_erg, color=color,
            linestyle='--', zorder=1)
    ax.plot(energies,
            bowtie['dnde_hi'] * e2 * MeV_to_erg, color=color,
            linestyle='--', zorder=1)

    ax.fill_between(energies,
                    bowtie['dnde_lo'] * e2 * MeV_to_erg,
                    bowtie['dnde_hi'] * e2 * MeV_to_erg,
                    alpha=0.5, color=color, zorder=-1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy [GeV]')
    ax.set_ylabel(r'$\nu f(\nu)$ [erg cm$^{-2}$ s$^{-1}$]')
    source = llh['config']['selection']['target']
    tmin = llh['config']['selection']['tmin']
    tmax = llh['config']['selection']['tmax']
    ax.text(0.2, 1.02, source,
            horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(0.5, 1.02, 'TS: {:.1f}'.format(llh['sources'][source]['ts']),
            horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 1.02, 'MJD: {:.1f} - {:.1f}'.format(MET_to_MJD(tmin), MET_to_MJD(tmax)),
            horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes)
    plt.grid(True)
    plt.tight_layout()
    fig.savefig(os.path.join(basepath, 'sed.pdf'),
                bbox_inches='tight')

def make_edges(data_f):
    header = data_f[0].header
    xmin = header['CRVAL1']-header['CRPIX1']*abs(header['CDELT1'])
    xmax = header['CRVAL1']+header['CRPIX1']*abs(header['CDELT1'])
    ymin = header['CRVAL2']-header['CRPIX2']*abs(header['CDELT2'])
    ymax = header['CRVAL2']+header['CRPIX2']*abs(header['CDELT2']) 
    xbins = header['NAXIS1']
    ybins = header['NAXIS2']
    return xmin, xmax, ymin, ymax, xbins, ybins


def make_ts_plot(basepath, srcs, mode='ts'):
    markers = ['o', 's', 'P', 'p', '*' , 'x', 'X', 'D', 4, 5, 6, 7]
    if mode == 'ts':
        inp = fits.open(os.path.join(basepath,'fit1_pointsource_powerlaw_2.00_tsmap.fits'))
    elif mode == 'resmap':
        inp = fits.open(os.path.join(basepath,'fit1_pointsource_powerlaw_2.00_residmap.fits'))
    srcs = np.load(srcs)[()]
    fig = plt.figure(figsize=figsize(0.4, 1.))
    plt.clf()
    ax=fig.add_axes((.0, .0,0.7,.7))
    xmin, xmax, ymin, ymax, xbins, ybins = make_edges(inp)
    X=np.linspace(xmin, xmax, xbins)[::-1]
    Y=np.linspace(ymin, ymax, ybins)
    xs, ys = np.meshgrid(X,Y)
    if mode == 'ts':
        Z=inp[2].data
        minmax = (0,9)
        cmap = ts_cmap
    elif mode == 'resmap':
        Z=inp[0].data
        minmax = (-3, 3)
        cmap = re_cmap
    cbar = ax.contourf(X,Y,Z,
                       levels=np.linspace(minmax[0], minmax[1], 500),
                       cmap=cmap) 
    levels=np.linspace(minmax[0], minmax[1], 4)
    CS = ax.contour(X,Y,Z, levels=levels,
                    colors='black', linewidths=(0.3,))
    ax.set_xlabel(r'R.A. (degrees)')
    ax.set_ylabel(r'Dec. (degrees)')
    plt.gca().invert_xaxis()
    ax.plot(inp[0].header['CRVAL1'], inp[0].header['CRVAL2'],
            marker='o', color='blue', ms=3, fillstyle='none')
    circle = plt.Circle((inp[0].header['CRVAL1'], inp[0].header['CRVAL2']), 1.0,
                        color='b', fill=False, linewidth=0.5)
    ax.add_artist(circle)
    for i, src in enumerate(zip(srcs['ra'], srcs['dec'])):
        ax.plot(src[0], src[1],
                marker=markers[i],
                linestyle = '',
                label=srcs['name'][i],
                color='k', ms=4)
    ax.set_xlim(inp[0].header['CRVAL1']+2.05, inp[0].header['CRVAL1']-2.05)
    ax.set_ylim(inp[0].header['CRVAL2']-2.05, inp[0].header['CRVAL2']+2.05)
    ax2=fig.add_axes((.0, .73,0.7,.05))
    plt_cbar = fig.colorbar(cbar, orientation="horizontal", cax=ax2,
                            ticks=np.arange(minmax[0], minmax[1] , 1))
    if mode == 'ts':
        plt_cbar.set_label(r'TS Value', labelpad=8)
    elif mode == 'resmap':
        plt_cbar.set_label(r'Significance [$\sigma$]', labelpad=8)
    plt_cbar.ax.xaxis.set_ticks_position('top')
    plt_cbar.ax.xaxis.set_label_position('top')
    ax2.tick_params(axis='x', which='major', pad=3)
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    plt.tight_layout()
    plt.savefig(os.path.join(basepath,'{}_map.pdf'.format(mode)),
                bbox_inches='tight')


def get_index(flux_dict):
    if flux_dict['SpectrumType']=='PowerLaw':
        return flux_dict['spectral_pars']['Index']['value']
    

def get_index_err(flux_dict):
    if flux_dict['SpectrumType']=='PowerLaw':
        return flux_dict['spectral_pars']['Index']['error']
    
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))
