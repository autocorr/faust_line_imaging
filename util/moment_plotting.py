#!/usr/bin/env python3
"""
Moment Plotting
===============
Create figure galleries of the moment maps for a FAUST source. This script is
not meant to be run under CASA but a user's Python v3 distribution. Plots can
be made using the :func:`plot_all_moments` function.
"""

import re
import warnings
from glob import glob
from copy import deepcopy
from pathlib import Path

import numpy as np
import scipy as sp
import pandas as pd
from skimage import morphology
from matplotlib import pyplot as plt
from matplotlib import (patheffects, colors)
from matplotlib.ticker import AutoMinorLocator

import aplpy
import radio_beam
from astropy.io import fits
from astropy import units as u
from astropy import (convolution, coordinates, wcs)


warnings.simplefilter('ignore')  # NOTE stop that matplotlib deprecation warning
plt.rc('font', size=10, family='serif')
plt.rc('text', usetex=True)
plt.rc('xtick', direction='out')
plt.rc('ytick', direction='out')


CLR_CMAP = plt.cm.Spectral_r
CLR_CMAP.set_bad('0.5', 1.0)
HOT_CMAP = plt.cm.afmhot
HOT_CMAP.set_bad('0.5', 1.0)  # 0.2 for darker
#RDB_CMAP = plt.cm.coolwarm_r  # also seismic_r
RDB_CMAP = plt.cm.RdBu_r
RDB_CMAP.set_bad('0.5', 1.0)
TAU_CMAP = plt.cm.YlGnBu_r
TAU_CMAP.set_bad('0.5', 1.0)
VIR_CMAP = plt.cm.viridis
VIR_CMAP.set_bad('0.5', 1.0)

VELOS = {
        'BHB07-11':          3.6,
        'CB68':              5.0,
        'Elias29':           5.0,
        'GSS30':             3.5,
        'IRAS_15398-3359':   5.3,
        'IRS63':             0.0,
        'L1527':             5.9,
        'L1551_IRS5':        6.5,
        'L483':              5.5,
        'NGC1333_IRAS4A1A2': 7.0,
        'NGC1333_IRAS4C':    7.6,
        'R_CrA_IRS7B':       6.2,
        'VLA1623A':          4.0,
}


def open_fits(filen, relative=True):
    path = Path(filen).expanduser()
    return fits.open(path)


def default_mask(hdul):
    data = deepcopy(hdul[0].data)
    mask = np.nan_to_num(data) > 0
    mask = morphology.remove_small_objects(mask, min_size=64)
    kernel = morphology.disk(1)
    mask = morphology.opening(mask.squeeze(), kernel)
    return mask


def get_map_center(hdul):
    shape = hdul[0].data.shape
    im_wcs = wcs.WCS(hdul[0].header)
    decim_deg = im_wcs.all_pix2world(shape[0]/2, shape[1]/2, 0)
    ra, dec = decim_deg[0].item(), decim_deg[1].item()
    co = coordinates.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    return co


def convert_jy_to_k(hdul):
    header = hdul[0].header
    beam = radio_beam.Beam.from_fits_header(header)
    name = hdul[0].fileinfo()['file'].name
    freq = re.search('\d+\.\d+GHz', name).group().rstrip('GHz')
    freq = float(freq) * u.GHz
    conv = (1 * u.Jy).to(u.K, u.brightness_temperature(beam, freq))
    hdul[0].data = hdul[0].data * conv.value
    return conv


def add_label(label, xy, fontsize=10):
    txt = plt.annotate(label, xy=xy, xycoords='axes fraction',
        fontsize=fontsize)
    txt.set_path_effects([patheffects.withStroke(linewidth=4.5,
        foreground='w')])
    return txt


def save_figure(filen, do_eps=False):
    exts = ['png', 'pdf']
    if do_eps:
        exts.append('eps')
    for ext in exts:
        path = 'plots' / Path('{0}.{1}'.format(filen, ext))
        plt.savefig(str(path), dpi=300)
        print('-- {0} saved'.format(ext))
    plt.close('all'); plt.cla(); plt.clf()


def fix_label_ticks(gc, co, radius=1.5/60):
    gc.recenter(co.fk5.ra.value, co.fk5.dec.value, radius=radius)
    gc.axis_labels.set_ypad(-7)
    gc.ticks.set_color('black')
    gc.tick_labels.set_xformat('hh:mm:ss')
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.axis_labels.set_font(size='small')
    gc.tick_labels.set_font(size='x-small')


def add_beam(gc, corner='top left', facecolor='black',
        edgecolor='none'):
    gc.add_beam()
    gc.beam.set_corner(corner)
    gc.beam.set_color(facecolor)
    gc.beam.set_edgecolor(edgecolor)
    gc.beam.set_linewidth(1.0)


def add_colorbar(gc, box, label=None, top=False):
    gc.add_colorbar(box=box, pad=0.2, width=0.1, axis_label_text=label)
    gc.colorbar.set_font(size='x-small')
    gc.colorbar._colorbar.ax.tick_params(direction='in')
    if top:
        gc.colorbar.set_location('top')
        gc.colorbar.set_width(0.1)
        gc.colorbar.set_pad(0.12)
        gc.colorbar._colorbar.ax.tick_params(direction='in')
        gc.colorbar._colorbar.ax.xaxis.set_ticks_position('top')


def add_primary_beam_contour(gc, hdul, power=0.205):
    gc.show_contour(hdul, dimensions=[0,1], slices=[100,0], levels=[power],
            linewidths=0.5, linestyles='dashed', colors='0.25',
            convention='calabretta')


def zero_out_corners(hdu):
    hdu[0].data[ 0, 0] = 0
    hdu[0].data[ 0,-1] = 0
    hdu[0].data[-1, 0] = 0
    hdu[0].data[-1,-1] = 0


def plot_four_moments(field, spw_label):
    """
    Parameters
    ----------
    field : str
        Target source name, e.g. 'CB68'.
    spw_label : str
        SPW label, e.g. '244.936GHz_CS'.
    """
    # FIXME due to issues with having 3 axes from STOKES/RA/DEC, `aplpy==2.0.1`
    # does not work and must currently fall-back to `aplpy==1.1.1` that does
    # not wrap the new WCSAxes module in astropy
    filen = 'moments/{0}_{1}'.format(field, spw_label)
    img0 = open_fits(filen + '.mom0_cmask.pbcor.fits')
    img1 = open_fits(filen + '.max.pbcor.fits')
    img2 = open_fits(filen + '.mom1.fits')
    img3 = open_fits(filen + '.mom2.fits')
    convert_jy_to_k(img0)
    convert_jy_to_k(img1)
    for hdul in (img0, img1, img2, img3):
        data = hdul[0].data
        mask = default_mask(hdul).reshape(data.shape)
        data[~mask] = np.nan
        zero_out_corners(hdul)
    for hdul in (img0, img1):
        data = hdul[0].data
        data[~np.isfinite(data)] = 0.0
    fig = plt.figure(figsize=(8, 6.5))
    split_size = 0.02
    x_margin_size = 0.15
    y_margin_size = 0.075
    delta = 0.06
    x0 = x_margin_size - delta
    y0 = y_margin_size
    dx = (1 - 2 * x_margin_size - split_size) / 2
    dy = (1 - 2 * y_margin_size - split_size) / 2
    x1 = x0 + dx + split_size + 1.5 * delta
    y1 = y0 + dy + split_size + 0.015
    labelx = 0.035
    labely = 0.035
    xlabel = r'$\alpha \ (\mathrm{J2000})$'
    ylabel = r'$\delta \ (\mathrm{J2000})$'
    # Moment 0, integrated intensity
    gc0 = aplpy.FITSFigure(img0, figure=fig, subplot=[x0, y1, dx, dy])
    gc0.show_colorscale(vmin=0, vmax=15, cmap=HOT_CMAP)
    gc0.tick_labels.hide_x()
    gc0.axis_labels.hide_x()
    gc0.axis_labels.set_ytext(ylabel)
    add_label(r'Moment 0 $(\int I \, dv)$', (labelx, labely))
    add_colorbar(gc0, [x0+dx+0.01, y1, 0.02, dy],
        label=r'$\int I \, dv \ \ [\mathrm{K \, km \, s^{-1}}]$')
    add_beam(gc0, facecolor='white')
    # Maximum, peak intensity
    gc1 = aplpy.FITSFigure(img1, figure=fig, subplot=[x1, y1, dx, dy])
    gc1.show_colorscale(vmin=0, vmax=10, cmap=HOT_CMAP, stretch='power',
            exponent=1)
    add_label(r'Maximum $(I_\mathrm{pk})$', (labelx, labely))
    add_colorbar(gc1, [x1+dx+0.01, y1, 0.02, dy],
        label=r'$T_\mathrm{pk} \ \ [\mathrm{K}]$')
    gc1.axis_labels.hide()
    gc1.tick_labels.hide()
    add_beam(gc1, facecolor='white')
    # Moment 1, intensity weighted velocity
    gc2 = aplpy.FITSFigure(img2, figure=fig, subplot=[x0, y0, dx, dy])
    vcen = VELOS[field]
    gc2.show_colorscale(vmin=vcen-1.5, vmax=vcen+1.5, cmap=RDB_CMAP)
    gc2.axis_labels.set_xtext(xlabel)
    gc2.axis_labels.set_ytext(ylabel)
    add_label(r'Moment 1 $(v_\mathrm{lsr})$', (labelx, labely))
    add_colorbar(gc2, [x0+dx+0.01, y0, 0.02, dy],
        label=r'$v_\mathrm{lsr} \ \ [\mathrm{km \, s^{-1}}]$')
    add_beam(gc2)
    # Moment 2, velocity dispersion
    gc3 = aplpy.FITSFigure(img3, figure=fig, subplot=[x1, y0, dx, dy])
    gc3.show_colorscale(vmin=0.15, vmax=0.8, cmap=CLR_CMAP)
    add_label(r'Moment 2 $(\sigma_v)$', (labelx, labely))
    add_colorbar(gc3, [x1+dx+0.01, y0, 0.02, dy],
        label=r'$\sigma_v \ \ [\mathrm{km \, s^{-1}}]$')
    gc3.tick_labels.hide_y()
    gc3.axis_labels.hide_y()
    gc3.axis_labels.set_xtext(xlabel)
    add_beam(gc3)
    # Common axis options
    center = get_map_center(img0)
    for gc in (gc0, gc1, gc2, gc3):
        fix_label_ticks(gc, center, radius=20*u.arcsec.to('deg'))
    out_filen = '{0}_{1}_moment_four_panel'.format(field, spw_label)
    save_figure(out_filen)


def plot_all_moments(field):
    """
    Parameters
    ----------
    field : str
        Target source name, e.g. 'CB68'.
    """
    for filen in glob(f'moments/{field}_*_clean.mom0.pbcor.fits'):
        basename = Path(filen).stem
        spw_label = re.search('.*clean', basename).group()[len(field)+1:]
        plot_four_moments(field, spw_label)


