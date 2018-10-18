# coding: utf-8

import sys
sys.path.append('/home/ga53lag/Software/python_scripts/')
from fancy_plot import *
import numpy as np
import os
from myfunctions import dict_to_nparray, MET_to_MJD, GreatCircleDistance
from astropy.time import Time
import pyfits as fits
import warnings
from matplotlib.colors import LinearSegmentedColormap
ts_cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', '#800000'])
re_cmap = LinearSegmentedColormap.from_list('mycmap2', ['#67a9cf', '#f7f7f7', '#ef8a62'])
MeV_to_erg = 1.60218e-6

def make_lc_plot(basepath):
    lc_dict = dict()
    keys = ['ts', 'dgamma','gamma',  'flux_ul95',  'flux_err', 'flux']
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
    print lc_arr.dtype.fields
    print lc_arr
    mask  = (np.abs(lc_arr['ts'])>4) & ((lc_arr['dgamma']/lc_arr['gamma'])<0.5)

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
    ax1.set_xlabel('Time (MJD)')

    ## Spectral Index Axis
    ax2=fig.add_axes((.0, .6,1.,.4))
    av_gamma, std_gamma = weighted_avg_and_std(lc_arr['gamma'],
                                               weights=1./lc_arr['dgamma'])
    ax2.axhline(av_gamma,linestyle='--', color='grey',
                alpha=0.85, zorder = -1, linewidth=0.9)
    print lc_arr['tmid'][mask]
    print lc_arr['gamma'][mask]
    ax2.errorbar(lc_arr['tmid'][mask], lc_arr['gamma'][mask],
                 yerr=lc_arr['dgamma'][mask],
                 xerr=lc_arr['bin_len'][mask]/2., linestyle='')
    ax2.set_ylabel('Index', labelpad=5)
    ax2.set_ylim(1.1,4)
    ax2.set_xticks([])
    ax2.set_xlim(ax1.get_xlim()[0], ax1.get_xlim()[1])
    ax2.text(0.8, 1.1, source,
        horizontalalignment='center',
        verticalalignment='center', transform=ax2.transAxes)
    plt.savefig(os.path.join(basepath, 'lightcurve.pdf'),
                bbox_inches='tight')
    return

def make_sed_plot(basepath):
    try:
        sed = np.load(os.path.join(basepath, 'sed.npy'))[()]
        bowtie = np.load(os.path.join(basepath, 'bowtie.npy'))[()]
        llh = np.load(os.path.join(basepath, 'llh.npy'))[()]
    except Exception as inst:
        warnings.warn('SED files not found')
        print(inst)
        return
    source = llh['config']['selection']['target']
    ts = llh['sources'][source]['ts']
    fig, ax = newfig(0.9)
    ax.set_xscale('log')
    ax.set_yscale('log')   
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

    if ts>4:
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
    ax.errorbar(x[~m], y[~m], xerr=xerr1,
                yerr=(yerr_lo[~m], yerr_hi[~m]),
                linestyle='', color='k', fmt='o', zorder=2)
    factor = np.log10(ax.get_ylim()[1]) - np.log10(ax.get_ylim()[0])
    serr = 10**np.log10(yul[m]) - 10**(np.log10(yul[m])-0.2/4.*factor)
    ax.errorbar(x[m], yul[m], xerr=xerr0,
                yerr=serr,  uplims=True,
                color='k', fmt='o', zorder=2)


    ax.set_xlabel('Energy [GeV]')
    ax.set_ylabel(r'$\nu f(\nu)$ [erg cm$^{-2}$ s$^{-1}$]')
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


def make_ts_plot(basepath, srcs, vou_cand, mode='tsmap'):
    markers = ['o', 's', 'P', 'p', '*' , 'x', 'X', 'D', 4, 5, 6, 7]
    ROI = 130./60.
    fname = 'fit1_pointsource_powerlaw_2.00_{}.fits'.format(mode)
    try:
        inp = fits.open(os.path.join(basepath, fname))
    except Exception as inst:
        warnings.warn('Files for {} not found'.format(mode))
        print(inst)
        return
    srcs = np.load(srcs)[()]
    cand_pos = np.genfromtxt(vou_cand)

    cand = {'ra': cand_pos[:,0], 'dec': cand_pos[:,1]}
    distances = np.ones(len(cand['ra']))
    for i in range(len(cand['ra'])):
        tdist= [GreatCircleDistance(
                  cand['ra'][i], cand['dec'][i], srcs['ra'][j],
                  srcs['dec'][j], unit='deg') for j in range(len(srcs['ra']))]
        distances[i] = np.degrees(np.min(tdist))
    inds =  np.argsort(distances)[-len(srcs['ra'])-1:]
    cand['ra'] = cand['ra'][inds]
    cand['dec'] = cand['dec'][inds]

    fig = plt.figure(figsize=figsize(0.4, 1.))
    plt.clf()
    ax=fig.add_axes((.0, .0,0.7,.7))
    xmin, xmax, ymin, ymax, xbins, ybins = make_edges(inp)
    X=np.linspace(xmin, xmax, xbins)[::-1]
    Y=np.linspace(ymin, ymax, ybins)
    xs, ys = np.meshgrid(X,Y)
    if mode == 'tsmap':
        Z=inp[2].data
        minmax = (0,25)
        ticks = 5
        cmap = ts_cmap
    elif mode == 'residmap':
        Z=inp[0].data
        minmax = (-5, 5)
        ticks = 2 
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
    circle = plt.Circle((inp[0].header['CRVAL1'], inp[0].header['CRVAL2']), 1.5,
                        color='b', fill=False, linewidth=0.5)
    ax.add_artist(circle)
    for i, src in enumerate(zip(srcs['ra'], srcs['dec'])):
        ax.plot(src[0], src[1],
                marker=markers[i],
                linestyle = '',
                label=srcs['name'][i],
                color='k', ms=4)
    ax.plot(cand['ra'], cand['dec'], color='k', ms=4, marker="8",
            linestyle = '', fillstyle='none', label='VOU Sources')
    ax.set_xlim(inp[0].header['CRVAL1']+ROI, inp[0].header['CRVAL1']-ROI)
    ax.set_ylim(inp[0].header['CRVAL2']-ROI, inp[0].header['CRVAL2']+ROI)
    ax2=fig.add_axes((.0, .73,0.7,.05))
    plt_cbar = fig.colorbar(cbar, orientation="horizontal", cax=ax2,
                            ticks=np.arange(minmax[0], minmax[1] , ticks))
    if mode == 'tsmap':
        plt_cbar.set_label(r'TS Value', labelpad=8)
    elif mode == 'residmap':
        plt_cbar.set_label(r'Significance [$\sigma$]', labelpad=8)
    plt_cbar.ax.xaxis.set_ticks_position('top')
    plt_cbar.ax.xaxis.set_label_position('top')
    ax2.tick_params(axis='x', which='major', pad=3)
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    plt.tight_layout()
    plt.savefig(os.path.join(basepath,'{}.png'.format(mode)),
                bbox_inches='tight', dpi=300)


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
