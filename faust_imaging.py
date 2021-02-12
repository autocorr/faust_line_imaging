#!/usr/bin/env python
"""
===================
ALMA FAUST Pipeline
===================
Line imaging pipeline for the ALMA FAUST Large Program. Please see the README
file for further documentation.

NOTE: To use this script within an interactive session of CASA (Python v2), run
with `execfile`.

Author:`Brian Svoboda`
Contributors:`Claire Chandler`
License:MIT
Year:2020
"""
from __future__ import (print_function, division)

import os
import shutil
import datetime
from glob import glob
from collections import (OrderedDict, Iterable)

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib import patheffects as path_effects

from cleanhelper import cleanhelper


# matplotlib configuration settings
plt.rc('font', size=10, family='serif')
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
plt.ioff()  # turn off interactive GUI

# Global paths
DATA_DIR = '/lustre/aoc/users/cchandle/FAUST/2018.1.01205.L/completed_SBs/'
PROD_DIR = '/lustre/aoc/users/bsvoboda/faust/faust_alma/data/'
IMAG_DIR = os.path.join(PROD_DIR, 'images/')
MOMA_DIR = os.path.join(PROD_DIR, 'moments/')
PLOT_DIR = os.path.join(PROD_DIR, 'plots/')

# Median-absolute-deviation conversion factor to Gaussian RMS
MAD_TO_RMS = 1 / (np.sqrt(2) * sp.special.erfinv(0.5))  # approx. 1.4826
# Large `niter` value, hopefully not reached in actually executions.
NITERMAX = int(1e7)
# Primary beam limit to image down to. Values in the image are mainly dominated
# by 7m data outside of a pblimit=0.07 and will have a non-uniformly sized
# synthesized beam.
PBLIMIT = 0.2
# Sampling factor to determine cell size relative to target resolution (~0.3as,
# see values `ALL_TARGETS`). A relatively large value of 10 pix per synthesized
# HPBW is used for consistency with the continuum images.
OVERSAMPLE_FACT = 10
# Line cut-out window in velocity (full width)
LINE_VWIN = '20km/s'
# Set of briggs uv-weightings to use by default in pipeline tasks
WEIGHTINGS = (0.5,)
# Default extension name for dirty images
DIRTY_EXT = 'dirty'
# Default extension name for unmasked clean model
NOMSK_EXT = 'nomask'
# perchanelweightdensity parameter in tclean, changing to False slightly
# improves the resolution at the expense of sensitivity and uniform RMS.
PERCHANWT = False
# See guide for discussion for of automasking parameters:
#   https://casaguides.nrao.edu/index.php/Automasking_Guide
AUTOM_KWARGS = {
        'usemask': 'auto-multithresh',
        'noisethreshold': 5.0,
        'sidelobethreshold': 2.0,
        'lownoisethreshold': 1.0,
        'minbeamfrac': 0.001,
        'negativethreshold': 5.0,
        'growiterations': 1000,
        'dogrowprune': False,
        'fastnoise': False,
        'verbose': True,
}


class Target(object):
    def __init__(self, name, res, vsys):
        self.name = name
        self.res = res
        self.vsys = vsys

    @property
    def vstart(self):
        half_vwin = qa.div(qa.quantity(LINE_VWIN), 2)
        vsys = qa.quantity(self.vsys, 'km/s')
        voffset = qa.sub(vsys, half_vwin)
        return qa.tos(qa.convert(voffset, 'km/s'))


ALL_TARGETS = { t.name: t for t in [
        Target('BHB07-11',          0.35, 3.6),
        Target('CB68',              0.36, 5.0),
        Target('Elias29',           0.36, 5.0),
        Target('GSS30',             0.36, 3.5),
        Target('IRAS_15398-3359',   0.32, 5.3),
        Target('IRS63',             0.34, 0.0),
        Target('L1527',             0.36, 5.9),
        Target('L1551_IRS5',        0.34, 6.5),
        Target('L483',              0.25, 5.5),
        Target('NGC1333_IRAS4A1A2', 0.21, 7.0),
        Target('NGC1333_IRAS4C',    0.21, 7.6),
        Target('R_CrA_IRS7B',       0.37, 6.2),
        Target('VLA1623A',          0.36, 4.0),
]}
ALL_FIELD_NAMES = ALL_TARGETS.keys()


class DataSet(object):
    ms_fmt = DATA_DIR + '{0}-Setup{1}/uid___*_target_lines_self_calibrated_continuum_subtracted_aligned.ms'
    low_freqs = {1: 217, 2: 245, 3: 93}  # GHz
    high_freqs = {1: 235, 2: 262, 3: 108}  # GHz

    def __init__(self, field, setup=1, kind='joint'):
        """
        Parameters
        ----------
        field : str
            FAUST target field name.
        setup : int
            FAUST Setup index number:
                1 : Band 6, 220 GHz
                2 : Band 6, 250 GHz
                3 : Band 3,  90 GHz
        kind : str
            Datset descriptor for which measurement set files to use. Valid
            values include: ('joint', '12m', '7m').
        """
        assert field in ALL_FIELD_NAMES
        assert setup in (1, 2, 3)
        kind = kind.lower()
        self.field = field
        self.target = ALL_TARGETS[field]
        self.setup = setup
        self.kind = kind
        vis_filen = glob(self.ms_fmt.format(field, setup))
        if kind == 'joint':
            self.vis = vis_filen
        elif kind == '12m':
            # use msmd to get array diameter size
            self.vis = [
                    name for name in vis_filen
                    if ms_contains_diameter(name, diameter=12)
            ]
        elif kind == '7m':
            self.vis = [
                    name for name in vis_filen
                    if ms_contains_diameter(name, diameter=7)
            ]
        else:
            raise ValueError('Invalid dataset descriptor: "{0}"'.format(kind))
        self.check_if_product_dirs_exist()

    @property
    def cell_12m(self):
        cell = self.target.res / OVERSAMPLE_FACT
        return '{0:.4f}arcsec'.format(cell)

    @property
    def cell_7m(self):
        hi_nu = self.high_freqs[self.setup]
        freq = qa.quantity(hi_nu, 'GHz')
        diam = '50m'  # max-baseline for ACA
        angle = self.calc_res(freq, diam)
        cell = angle['value'] / OVERSAMPLE_FACT
        return '{0:.4f}arcsec'.format(cell)

    @property
    def cell(self):
        if self.kind in ('joint', '12m'):
            return self.cell_12m
        elif self.kind == '7m':
            return self.cell_7m
        else:
            raise ValueError

    @property
    def gridder(self):
        return 'mosaic' if self.kind == 'joint' else 'standard'

    @property
    def pblimit(self):
        return PBLIMIT if self.kind == 'joint' else -0.01

    @property
    def imsize(self):
        if self.kind in ('joint', '12m'):
            ant_diam = '12m'
        elif self.kind == '7m':
            ant_diam = '7m'
        else:
            raise ValueError
        diam = qa.div(ant_diam, 1.13)  # account for tapered illumination
        freq = qa.quantity(self.low_freqs[self.setup], 'GHz')
        angle = self.calc_res(freq, diam)
        angle = qa.mul(angle, 2.00)  # 1.1 * [full-width at 10% maximum]
        pix_width = qa.convert(qa.div(angle, self.cell), '')['value']
        eff_width = cleanhelper.getOptimumSize(int(pix_width))
        return [eff_width, eff_width]

    def check_if_product_dirs_exist(self):
        field_dir = os.path.join(IMAG_DIR, self.field)
        try:
            os.makedirs(field_dir)
        except OSError:
            if not os.path.isdir(field_dir):
                raise

    def calc_res(self, freq, diam):
        """Note that no factor of 1.028 is applied for FWHM of Airy pattern, just l/D."""
        assert qa.isquantity(freq)
        assert qa.isquantity(diam)
        wavel = qa.convertfreq(freq, 'm')
        angle = qa.convert(qa.mul(qa.div(wavel, diam), '1rad'), 'arcsec')
        assert qa.isangle(angle)
        return angle


def ms_contains_diameter(ms_filen, diameter=None, epsilon=0.5):
    """
    Determine whether an MS contains an antenna with a given diameter.

    Parameters
    ----------
    diameter : number
        Antenna diameter in meters.
    epsilon : number, default 0.5
        Fudge factor for querying the antennas min/max diameter.
    """
    mindiameter = '{0}m'.format(diameter-epsilon)
    maxdiameter = '{0}m'.format(diameter+epsilon)
    msmd.open(ms_filen)
    ant_ids = msmd.antennaids(mindiameter=mindiameter, maxdiameter=maxdiameter)
    msmd.close()
    # NOTE this will be incorrect if the MS files have been concat
    return ant_ids.size > 0


class Spw(object):
    line_vwin = '20km/s'

    def __init__(self, setup, restfreq, mol_name, name, ms_restfreq, spw_id, ot_name, nchan, chan_width,
            tot_bw):
        """
        Parameters
        ----------
        setup : int
            Setup ID number (1, 2, 3)
        restfreq : str
            Rest frequency of primary targeted line in the SPW, e.g. "93.17GHz".
        mol_name : str
            Molecule name of the primary targeted line in the SPW.
        name : str
            Name of the spectral line.
        ms_restfreq : str
            Line rest frequency listed in the measurement set, e.g. "93.17GHz".
            Many of these values are not rest frequencies for specific molecular
            transitions but shifted values meant to center the bandpass.
        spw_id : int
            ID number of the spectral window
        ot_name : str
            Ending of the OT SPW label, e.g. "#BB_1#SW-01". These are useful
            for selecting the spectral window ID number if the correlator ID
            numbers are inconsistent across 7M/12M datasets.
        nchan : int
            Total number of channels (before flagging).
        chan_width : number
            Channel width in kHz (TOPO).
        tot_bw : number
            Total bandwidth in MHz (TOPO).
        """
        assert setup in (1, 2, 3)
        self.setup = setup
        self.restfreq = restfreq
        self.mol_name = mol_name
        self.name = name
        self.ms_restfreq = ms_restfreq
        self.spw_id = spw_id
        self.ot_name = ot_name
        self.nchan = nchan
        self.chan_width = chan_width
        self.tot_bw = tot_bw

    @property
    def win_chan(self):
        chan_fwidth = qa.quantity(self.chan_width, 'Hz')
        chan_vwidth = qa.convertdop(qa.div(chan_fwidth, self.restfreq), 'km/s')
        nchan = qa.convert(qa.div(LINE_VWIN, chan_vwidth), '')['value']
        return int(np.ceil(nchan))

    @property
    def short_restfreq(self):
        unit = 'GHz'
        freq = qa.convert(self.restfreq, unit)
        return '{0:.3f}{1}'.format(freq['value'], unit)

    @property
    def label(self):
        return '{0}_{1}'.format(self.short_restfreq, self.mol_name)


def parse_spw_info_from_ms(ms_filen, field='CB68', out_filen=None):
    # NOTE This function was only needed to write out the values that can now
    # be found in the `Spw` objects. The resource configuration should be the
    # same for all targets.
    msmd.open(ms_filen)
    # get source ID number
    field_id = msmd.fieldsforname(field)[0]
    # SPW data descriptor IDs just for science observations
    sci_dds = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
    # list of SPW names from the science SPW DD IDs
    spw_names = msmd.namesforspws(sci_dds)
    spw_names = ['#'.join(n.split('#')[2:]) for n in spw_names]
    spw_metadata = []
    for ii, dd in enumerate(sci_dds):
        # transitions
        trans = msmd.transitions(field_id, dd)[0]
        # restfreqs
        freq_q = msmd.restfreqs(field_id, dd)['0']['m0']
        freq_s = '{0}{1}'.format(freq_q['value'], freq_q['unit'])
        # nchan
        nchan = msmd.nchan(dd)
        # channel resolution
        chanres = abs(msmd.chanres(dd).mean())  # in Hz
        # bandwidth
        bandwidth = msmd.bandwidths(dd)  # in Hz
        spw_metadata.append(
                (trans, freq_s, dd, spw_names[ii], nchan, chanres, bandwidth)
        )
    msmd.close()
    if out_filen is not None:
        with open(out_filen, 'w') as f:
            for meta in spw_metadata:
                f.write("Spw{0},\n".format(str(meta)))
    return spw_metadata


def write_all_spw_info():
    for setup_ix in range(1, 4):
        log_post('-- Writing metadata for Setup-{0}'.format(setup_ix))
        dset = DataSet(setup=setup_ix)
        filen = dset.vis[-1]
        out_filen = 'spw_info_Setup{0}.txt'.format(setup_ix)
        parse_spw_info_from_ms(filen, out_filen=out_filen)


def spw_list_to_dict(spws):
    return OrderedDict((s.label, s) for s in spws)


# NOTE These values are copy-pasted from the files generated by `write_all_spw_info`.
#      They should be equivalent for all FAUST targets. The SPW ID numbers are for
#      the 12m EBs, but these are only for reference. Actual SPW ID's are retrieved
#      from each MS based on the "BB" name.
# SPWs for Setup 1 (Band 6 "a", 220 GHz)
SPW_S1 = spw_list_to_dict([
    Spw(1, '216.112580GHz', 'DCOp', 'DCO__v_0_J_3_2__HCOOCH3_19_2_18__18_2_17__A(ID=0)', '2.16112582e+11Hz', 25, 'BB_1#SW-01#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '216.562710GHz', 'NH2D', 'NH2D_3_2_2_0s_3_1_2_0a__CH3CHO_v_t_0_11_1_10__10_1_9___E(ID=0)', '2.16562716e+11Hz', 27, 'BB_1#SW-02#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '217.104980GHz', 'SiO', 'SiO_v_0_5_4(ID=0)', '2.17104919e+11Hz', 29, 'BB_1#SW-03#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '217.822150GHz', 'c-C3H2', 'c_C3H2_v_0_6_0_6__5_1_5___6_1_6__5_0_5_(ID=0)', '2.17822148e+11Hz', 31, 'BB_1#SW-04#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '218.222190GHz', 'H2CO', 'H2CO_3_0_3__2_0_2_(ID=0)', '2.18222192e+11Hz', 33, 'BB_2#SW-01#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '218.440050GHz', 'CH3OH', 'CH3OH_v_t_0_4_2_2__3_1_2___NH2CHO_10_1_9__9_1_8___H2CO_3_2_2__2_2_1_(ID=0)', '2.18440063e+11Hz', 35, 'BB_2#SW-02#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '219.560358GHz', 'C18O', 'C18O_2_1(ID=0)', '2.19560354e+11Hz', 37, 'BB_2#SW-03#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '219.949442GHz', 'SO', 'SO_3__v_0_6_5__5_4_(ID=0)', '2.19949442e+11Hz', 39, 'BB_2#SW-04#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '231.060983GHz', 'OCS', 'OCS_19_18(ID=0)', '2.31060993e+11Hz', 41, 'BB_3#SW-01#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '231.220686GHz', '13CS', '13CS_v_0_5_4(ID=0)', '2.31220685e+11Hz', 43, 'BB_3#SW-02#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '231.321828GHz', 'N2Dp', 'DN2__J_3_2__CH3CHO_4_lines(ID=0)', '2.313218283e+11Hz', 45, 'BB_3#SW-03#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '231.410234GHz', 'D2CO', 'D2CO_4_0_4__3_0_3___CH3CHO_12_5_8__11_5_7__E(ID=0)', '2.31410234e+11Hz', 47, 'BB_3#SW-04#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(1, '233.795750GHz', 'cont', 'continuum(ID=0)', '2.3379575e+11Hz', 49, 'BB_4#SW-01#FULL_RES', 3840, 976562.5, 1875000000.0),
])

# SPWs for Setup 2 (Band 6 "b", 260 GHz)
SPW_S2 = spw_list_to_dict([
    Spw(2, '243.915826GHz', 'CH3OH', 'CH3OH_v_t_0_5_1_4__4_1_3____(ID=0)', '2.43915788e+11Hz', 25, 'BB_1#SW-01#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '244.048592GHz', 'H2CS', 'H2CS_7_1_6__6_1_5_(ID=0)', '2.44048504e+11Hz', 27, 'BB_1#SW-02#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '244.935560GHz', 'CS', 'CS_v_0_5_4(ID=0)', '2.44935557e+11Hz', 29, 'BB_1#SW-03#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '245.606320GHz', 'HC3N', 'HC3N_v_0_J_27_26__SO2_10_3_7__10_2_8_(ID=0)', '2.4560632e+11Hz', 31, 'BB_1#SW-04#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '246.700000GHz', 'cont', 'continuum(ID=0)', '2.467e+11Hz', 33, 'BB_2#SW-01#FULL_RES', 1920, 1128906.25, 1875000000.0),
    Spw(2, '260.189080GHz', 'NH2CHO', 'NH2CHO_12_2_10__11_2_9_(ID=0)', '2.60189848e+11Hz', 35, 'BB_3#SW-01#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '260.255080GHz', 'H13COp', 'HCOOCH3_21_3_18__20_3_17_E__A__H13CO__J_3_2(ID=0)', '2.6025508e+11Hz', 37, 'BB_3#SW-02#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '261.687366GHz', 'CH2DOH', 'CH2DOH_5_2_4__5_1_5___e0__HCOOCH3_21_7_14__20_7_13__E(ID=0)', '2.61692e+11Hz', 39, 'BB_3#SW-03#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '262.004260GHz', 'CCH', 'CCH_v_0_N_3_2__J_7_2_5_2__F_4_3__3_2(ID=0)', '2.6200426e+11Hz', 41, 'BB_3#SW-04#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '257.895673GHz', 'CH2DOH', 'CH2DOH_4_2_3__3_1_3___e1(ID=0)', '2.578956727e+11Hz', 43, 'BB_4#SW-01#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '258.255826GHz', 'SO', 'SO_3__v_0_6_6__5_5_(ID=0)', '2.582558259e+11Hz', 45, 'BB_4#SW-02#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '258.548819GHz', 'CH3OCH3', 'CH3OCH3_14_1_14__13_0_13__AE__EA__EE__AA(ID=0)', '2.58548819e+11Hz', 47, 'BB_4#SW-03#FULL_RES', 480, 141113.28125, 58593750.0),
    Spw(2, '259.034910GHz', 'HDCO', 'HDCO_4_2_2__3_2_1_(ID=0)', '2.5903491e+11Hz', 49, 'BB_4#SW-04#FULL_RES', 480, 141113.28125, 58593750.0),
])

# SPWs for Setup 3 (Band 3 "a", 90 GHz)
SPW_S3 = spw_list_to_dict([
    Spw(3,  '93.180900GHz', 'N2Hp', 'N2H__v_0_J_1_0(ID=3925982)', '93180900000.0Hz', 25, 'BB_1#SW-01#FULL_RES', 960, 70556.640625, 58593750.0),
    Spw(3,  '94.405163GHz', '13CH3OH', '13CH3OH_v_t_0_2__1_2__1__1_1___2_0_2__1_0_1______2__0_2__1__0_1___2__1_1__1__1_0_(ID=0)', '94407129000.0Hz', 27, 'BB_1#SW-02#FULL_RES', 960, 70556.640625, 58593750.0),
    Spw(3,  '94.999908GHz', 'cont', 'Continuum(ID=0)', '94999908103.2Hz', 29, 'BB_2#SW-01#FULL_RES', 3840, 976562.5, 1875000000.0),
    Spw(3, '107.013803GHz', 'CH3OH', 'CH3OH_v_t_0_3_1_3__4_0_4____(ID=0)', '1.07013831e+11Hz', 31, 'BB_3#SW-01#FULL_RES', 960, 70556.640625, 58593750.0),
    Spw(3, '108.039986GHz', 'l-C3D', 'C3D_11_2_9_2_f__c_C3HD_5_5_1__5_4_2_(ID=0)', '1.08064e+11Hz', 33, 'BB_3#SW-02#FULL_RES', 960, 70556.640625, 58593750.0),
    Spw(3, '104.239295GHz', 'SO2', 'SO2_v_0_10_1_9__10_0_10_(ID=96928)', '1.042392952e+11Hz', 35, 'BB_4#SW-01#FULL_RES', 960, 70556.640625, 58593750.0),
    Spw(3, '105.799113GHz', 'H13CCCN', 'H13CCCN_J_12_11(ID=455557)', '1.05799113e+11Hz', 37, 'BB_4#SW-02#FULL_RES', 960, 70556.640625, 58593750.0),
])

ALL_SETUPS = (1, 2, 3)
SPWS_BY_SETUP = {1: SPW_S1, 2: SPW_S2, 3: SPW_S3}
ALL_SPWS = {}
ALL_SPWS.update(SPW_S1)
ALL_SPWS.update(SPW_S2)
ALL_SPWS.update(SPW_S3)
ALL_SPW_LABELS = ALL_SPWS.keys()


###############################################################################
# General utility functions
###############################################################################

def log_post(msg):
    """
    Post a message to the CASA logger, logfile, and stdout/console.
    """
    print(msg)
    casalog.post(msg, 'INFO', 'faust_pipe')


def check_delete_image_files(imagename, parallel=False, preserve_mask=False):
    """
    Check for and remove (if they exist) files created by clean such as '.flux',
    '.image', etc.
    NOTE this function has issues with deleting tables made by clean in
    parallel mode, probably because the directory structure is different.

    Parameters
    ----------
    imagename : str
        The relative path name to the files to delete.
    parallel : bool, default False
        rmtables can't remove casa-images made with parallel=True, they must be
        manually removed.
    preserve_mask : bool, default False
        Whether to preserve the `.mask` file extension
    """
    log_post(':: Check for and remove existing files')
    exts = [
        '.flux', '.pb', '.image', '.alpha', '.alpha.error', '.weight',
        '.model', '.pbcor', '.psf', '.sumwt', '.residual', '.flux.pbcoverage',
        '.pb.tt0', '.image.tt0', '.weight.tt0', '.model.tt0', '.psf.tt0',
        '.sumwt.tt0', '.residual.tt0',
    ]
    if not preserve_mask:
        exts += ['.mask']
    # CASA image table directories
    for ext in exts:
        filen = imagename + ext
        if os.path.exists(filen):
            if parallel:
                log_post('-- Hard delete {0}'.format(ext))
                shutil.rmtree(filen)
            else:
                log_post('-- Removing {0}'.format(filen))
                rmtables(filen)
    # "Cannot delete X because it's not a table" -> so hard delete
    for ext in ('.residual', '.workdirectory'):
        filen = imagename + ext
        if os.path.exists(filen):
            log_post('-- Hard delete {0}'.format(ext))
            shutil.rmtree(filen)


def delete_all_extensions(imagename, keep_exts=None):
    """
    Parameters
    ----------
    imagename : str
    keep_exts : None, iterable
        A list of extensions to keep, example: ['mask', 'psf']
    """
    for filen in glob(imagename+'.*'):
        if keep_exts is not None and any(filen.endswith(ext) for ext in keep_exts):
            continue
        if_exists_remove(filen)


def if_exists_remove(filen):
    if os.path.exists(filen):
        log_post(':: Removing {0}'.format(filen))
        try:
            rmtables(filen)
            # the file will still exist if rmtables failed
            shutil.rmtree(filen)
            log_post('-- Hard delete: {0}'.format(filen))
        except OSError:
            pass


def delete_workdir(imagename):
    workdir = '{0}.workdirectory'.format(imagename)
    if os.path.exists(workdir):
        shutil.rmtree(workdir)


def replace_existing_file_with_new(old_filen, new_filen):
    """
    Replace an existing file with a new or temporary one.
    `old_filen` will be removed and replaced by `new_filen`.
    """
    if not os.path.exists(new_filen):
        raise IOError('File does not exist: "{0}"'.format(new_filen))
    rmtables(old_filen)
    # If a parallel/virtual image, then will not be recognized as a table in
    # the call to `rmtables`
    if os.path.exists(old_filen):
        log_post('-- Hard delete! "{0}"'.format(old_filen))
        shutil.rmtree(old_filen)
    os.rename(new_filen, old_filen)


def export_fits(imagename, velocity=False, overwrite=True):
    log_post(':: Exporting fits')
    exportfits(imagename, imagename+'.fits', dropstokes=True,
            velocity=velocity, overwrite=overwrite)


def concat_parallel_image(imagename):
    """
    Create a contiguous image cube file from a 'virtual image' generated from a
    parallel `tclean` run. Generates a new image appended with '.concat'.
    """
    _, ext = os.path.splitext(imagename)
    ext = ext.lstrip('.')
    outfile = imagename + '.CONCAT_TEMP'
    # parallel images generated by tclean have names formatted as:
    #   name.ext/name.n1.ext
    infiles = glob(imagename+'/*.{0}'.format(ext))
    if len(infiles) == 0:
        log_post(':: Not a parallel image! Passing.')
        return
    log_post(':: Concatenating file "{0}"'.format(imagename))
    if_exists_remove(outfile)
    ia.open(imagename)
    im_tool = ia.imageconcat(
            outfile=outfile,
            infiles=infiles,
            reorder=True,
    )
    im_tool.close()
    ia.done()
    replace_existing_file_with_new(imagename, outfile)


def concat_parallel_all_extensions(imagebase):
    exts = ['image', 'mask', 'model', 'pb', 'psf', 'residual', 'sumwt', 'weight']
    for ext in exts:
        imagename = '{0}.{1}'.format(imagebase, ext)
        if os.path.exists(imagename):
            concat_parallel_image(imagename)


def calc_common_coverage_range(imagebase):
    """
    Calculate the common frequency coverage amongst a set of MSs/EBs. The
    `.sumwt` image extension is used from an existing "dirty" cube.

    Parameters
    ----------
    imagebase : str

    Returns
    -------
    start : str
        The start frequency in the LSRK frame.
    nchan : int
        The number of channels from `start` at the native spectral resolution.
    """
    sumwt_filen = '{0}.sumwt'.format(imagebase)
    # get sum-of-the-weights values
    ia.open(sumwt_filen)
    sumwt = ia.getchunk().squeeze()
    csys = ia.coordsys().torecord()
    ia.close()
    # determine lowest and highest channels to use
    med_sumwt = np.nanmedian(sumwt)
    good_chans = np.argwhere(sumwt > 0.95 * med_sumwt)
    ix_lo = good_chans.min()
    ix_hi = good_chans.max()
    nchan = abs(ix_hi - ix_lo)
    # convert indices to frequencies
    spec_sys = csys['spectral2']
    sunit = spec_sys['unit']
    crpix = spec_sys['wcs']['crpix']
    crval = spec_sys['wcs']['crval']
    cdelt = spec_sys['wcs']['cdelt']
    faxis = (np.arange(len(sumwt)) + crpix) * cdelt + crval
    start_val = faxis[ix_lo:ix_hi].min()
    start = '{0}{1}'.format(start_val, sunit)
    log_post('-- start frequency: {0}'.format(start))
    log_post('-- nchan: {0}'.format(nchan))
    return start, nchan


def primary_beam_correct(imagebase):
    log_post(':: Primary beam correcting image: {0}'.format(imagebase))
    impbcor(
            imagename=imagebase+'.image',
            pbimage=imagebase+'.pb',
            outfile=imagebase+'.image.pbcor',
            overwrite=True,
    )


def smooth_cube_to_common_beam(imagename):
    """
    Use the `imsmooth` task to convert a cube with per-plane beams into one
    with a single common beam.
    """
    log_post(':: Smoothing cube to common beam scale')
    outfile = imagename + '.common'
    imsmooth(imagename=imagename, kernel='commonbeam', outfile=outfile,
            overwrite=True)


def copy_pb_mask(imagename, pbimage):
    """
    Copy primary beam T/F mask back into the image after being removed by a
    task such as `imsmooth`.
    """
    log_post(':: Copying primary beam mask into: {0}'.format(imagename))
    makemask(mode='copy', inpimage=pbimage, inpmask=pbimage+':mask0',
            output=imagename+':mask0', overwrite=True)


def calc_rms_from_image(imagename, chan_start=None, chan_end=None):
    """
    Calculate RMS values from the scaled MAD of all channels (unless given a
    range) for the given full imagename.

    Parameters
    ----------
    imagename : str
        CASA Image name, e.g., "244.936GHz_CS_joint_0.5_dirty.image"
    chan_start : (int, None)
        Start channel. If None, use full channel range.
    chan_end : (int, None)
        End channel. If None, use full channel range.

    Returns
    -------
    rms : number
    """
    if chan_start is None or chan_end is None:
        chans = None
    else:
        assert chan_start < chan_end
        chans = '{0}~{1}'.format(chan_start, chan_end)
    # check if directory exists
    if not os.path.exists(imagename):
        raise IOError('File or directory not found: "{0}"'.format(imagename))
    # run imstat over the cube and extract the calculated RMS
    stats = imstat(imagename=imagename, axes=[0, 1], chans=chans)
    mad_arr = stats['medabsdevmed']
    rms = MAD_TO_RMS * np.nanmedian(mad_arr)
    return rms


def create_mask_from_threshold(infile, outfile, sigma):
    rms = calc_rms_from_image(infile)
    thresh = sigma * rms
    immath(imagename=infile, mode='evalexpr', outfile=outfile,
            expr='iif(IM0>{0},1,0)'.format(thresh))


def make_multiscale_joint_mask(imagename, sigma=5.0, mask_ang_scales=(0, 1, 3),
        overwrite=True):
    """
    Create a joint mask from multiple smoothed copies of an image. A file
    with the `.summask` extension is created from the union of an RMS threshold
    applied to each image. The RMS is re-computed for each smoothed cube.

    Parameters
    ----------
    imagename : str
        Image name base without the extension.
    sigma : number
        Multiple of the RMS to threshold each image.
    mask_ang_scales : Iterable(number)
        Gaussian kernel FWHM in units of arcseconds to convolve each image
        with (note: not in pixel units).
    overwrite : bool
        Whether to overwrite the smoothed mask files if they exist.

    Results
    -------
        Files are written for each smoothing scale:
            smoothed image: '<IMG>_smooth{.3f}.image' (excluding scale=0)
            masked image:   '<IMG>_smooth{.3f}.mask'
        and the joint or unioned mask file across scales is written to:
            joint mask:     '<IMG>.summask'
    """
    log_post(':: Creating threshold based mask')
    n_scales = len(mask_ang_scales)
    original_image = '{0}.image'.format(imagename)
    pb_image = '{0}.pb'.format(imagename)
    base_names = [
            '{0}_smooth{1:.3f}'.format(imagename, s)
            for s in mask_ang_scales
    ]
    mask_files = [s+'.mask' for s in base_names]
    image_files = [s+'.image' for s in base_names]
    for scale, smooth_image, mask_image in zip(mask_ang_scales, image_files, mask_files):
        assert scale >= 0
        if scale == 0:
            # For a smoothing scale of zero, use the original image.
            create_mask_from_threshold(original_image, mask_image, sigma)
            continue
        if overwrite:
            if_exists_remove(smooth_image)
            if_exists_remove(mask_image)
        # Smooth the image with a Gaussian kernel of the given scale.
        imsmooth(imagename=original_image, kernel='gauss',
                major='{0}arcsec'.format(scale),
                minor='{0}arcsec'.format(scale), pa='0deg', targetres=True,
                outfile=smooth_image)
        # PB mask is erased after smoothing, so copy it back in from the PB
        copy_pb_mask(smooth_image, pb_image)
        # Threshold the image to create a mask
        create_mask_from_threshold(smooth_image, mask_image, sigma)
    # Add all of the masks for an effective chained "binary OR" operation
    sum_expr = '+'.join([
            'IM{0}'.format(i)
            for i in range(n_scales)
    ])
    outfile = '{0}.summask'.format(imagename)
    if_exists_remove(outfile)
    immath(imagename=mask_files, mode='evalexpr', outfile=outfile,
            expr='iif({0}>0,1,0)'.format(sum_expr))
    log_post('-- joint mask file written to: {0}'.format(outfile))


def format_cube_robust(weighting):
    if isinstance(weighting, (int, float)):
        return weighting, 'briggs'
    else:
        # e.g., string type such as "natural"
        return None, weighting


def format_rms(rms, sigma=1, unit='Jy'):
    return '{0}{1}'.format(sigma*rms, unit)


def check_max_residual(imagebase, sigma=5.5):
    """
    Check the residual image cube for channels with deviations higher than a
    multiple of the global RMS.

    Parameters
    ----------
    imagebase : str
        Image name without extension.
    sigma : number
        Threshold to apply as a multiple of the RMS.
    """
    imagename = '{0}.residual'.format(imagebase)
    if not os.path.exists(imagename):
        raise OSError('File not found: {0}'.format(imagename))
    stats = imstat(imagename=imagename)  # over all channels
    rms = MAD_TO_RMS * stats['medabsdevmed'][0]
    v_max = stats['max'][0]
    v_min = stats['min'][0]
    max_resid = max(v_max, abs(v_min)) / rms
    if max_resid > sigma:
        log_post('!! WARNING: maximum residual threshold exceeded')
        log_post('-- maximum residual: {0}'.format(max_resid))
        log_post('-- residual threshold: {0} ({1}-sigma)'.format(sigma*rms, sigma))
        log_post('-- max: {0}, index: {1}'.format(v_max, stats['maxpos']))
        log_post('-- min: {0}, index: {1}'.format(v_min, stats['minpos']))


###############################################################################
# Line imaging
###############################################################################

class ImageConfig(object):
    scales = [0, 15, 45, 135]  # pixels; point, 1.5, 4.5 beam hpbw's (10 pix)
    mask_ang_scales = [0, 1, 3]  # arcsec
    smallscalebias = -1.0
    gain = 0.05
    cyclefactor = 2.0
    parallel = MPIEnvironment().is_mpi_enabled
    autom_kwargs = AUTOM_KWARGS.copy()

    def __init__(self, dset, spw, fullcube=True, weighting=0.5):
        """
        Configuration object for `tclean` related custom tasks. Parameters
        specify deconvolution and image-cube properties such as the weighting.

        Parameters
        ----------
        dset : DataSet
        spw : Spw
        fullcube : bool
            Image the full spectral window or a small window around the SPWs
            defined rest frequency.
        weighting : (str, number)
            The name of the weighting method (string) or briggs factor for briggs
            weighting (number)
        """
        assert dset.setup == spw.setup
        self.dset = dset
        self.spw = spw
        self.fullcube = fullcube
        # weighting method, must be valid for tclean
        robust, weighting_kind = format_cube_robust(weighting)
        self.weighting = weighting
        self.robust = robust
        self.weighting_kind = weighting_kind

    @classmethod
    def from_name(cls, field, label, kind='joint', **kwargs):
        spw = ALL_SPWS[label]
        dset = DataSet(field, setup=spw.setup, kind=kind)
        return cls(dset, spw, **kwargs)

    @property
    def dirty_imagebase(self):
        return self.get_imagebase(ext=DIRTY_EXT)

    @property
    def nomask_imagebase(self):
        return self.get_imagebase(ext=NOMSK_EXT)

    @property
    def rms(self):
        imagename = '{0}.image'.format(self.dirty_imagebase)
        rms = calc_rms_from_image(imagename)
        log_post('-- RMS {0} Jy'.format(rms))
        return rms

    @property
    def selected_start_nchan(self):
        if self.fullcube:
            # Restrict coverage to common coverage cross EBs.
            start, nchan = calc_common_coverage_range(self.dirty_imagebase)
        else:
            # Restrict coverage to range around source velocity.
            # NOTE This may produce undesirable behaviour if the line rest
            #      frequency is next to the edge of the SPW.
            start = self.dset.target.vstart
            nchan = self.spw.win_chan
        return start, nchan

    @property
    def spw_ids(self):
        """
        Determine spectral window ID numbers based on the (unique) Data
        Descriptor name found in the correlator configuration of the MS.
        This avoids indexing issues between mixed sets of 12m & 7m EBs.
        """
        ot_name = self.spw.ot_name
        spw_ids_for_vis = []
        for filen in self.dset.vis:
            msmd.open(filen)
            # SPW data descriptor IDs just for science observations
            sci_dds = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
            # list of SPW names from the science SPW DD IDs
            spw_names = msmd.namesforspws(sci_dds)
            msmd.close()
            # select IDs that have corresponding matching names
            matches = [
                    dd for dd, name in zip(sci_dds, spw_names)
                    if name.endswith(ot_name)
            ]
            if len(matches) == 0:
                raise ValueError(
                        'OT_NAME "{0}" not found in MS "{1}"'.format(ot_name, filen))
            elif len(matches) >= 2:
                raise ValueError(
                        'Multiple matches found, is MS "{0}" concatenated?'.format(filen))
            else:
                spw_ids_for_vis.extend(matches)
        return [str(n) for n in spw_ids_for_vis]

    def get_imagebase(self, ext=None):
        """
        Get the relative path name for an image basename following:
            "images/<FIELD>/<LABEL>_<ARRAY>_<WEIGHTING>[_<EXT>]"
        Note that this does not include image extensions such as ".image"

        Parameters
        ----------
        ext : (str, Iterable, object)
            If `ext` is None, then compose path of setup, spw, array.
            If `ext` is an iterable, append each by underscores.
            If other object, must have a `__repr__` method for representation string.
        """
        stem = 'images/{0}/{0}_{1}_{2}_{3}'.format(
                self.dset.field, self.spw.label, self.dset.kind, self.weighting)
        if ext is None:
            return stem
        elif isinstance(ext, Iterable) and not isinstance(ext, str):
            ext_items = '_'.join([str(i) for i in ext])
            return '{0}_{1}'.format(stem, ext_items)
        else:
            return '{0}_{1}'.format(stem, ext)

    def make_dirty_cube(self):
        dset, spw = self.dset, self.spw
        # Book-keeping before run
        imagename = self.dirty_imagebase
        log_post(':: Creating dirty cube ({0})'.format(imagename))
        # To create dirty cube simply perform zero iterations.
        niter = 0
        start = None
        nchan = -1
        # run tclean
        delete_all_extensions(imagename)
        tclean(
            vis=dset.vis,
            imagename=imagename,
            field=dset.field,
            spw=self.spw_ids,
            specmode='cube',
            chanchunks=-1,
            outframe='lsrk',
            veltype='radio',
            restfreq=spw.restfreq,
            nchan=nchan,
            start=start,
            imsize=dset.imsize,
            cell=dset.cell,
            # gridder parameters
            gridder=dset.gridder,
            weighting=self.weighting_kind,
            robust=self.robust,
            perchanweightdensity=PERCHANWT,
            # deconvolver parameters
            pblimit=dset.pblimit,
            niter=niter,
            interactive=False,
            parallel=self.parallel,
        )
        self.cleanup(imagename)

    def clean_line_nomask(self, sigma=4.5, scale_upper_limit=60):
        """
        Parameters
        ----------
        sigma : number
            Global threshold to clean down to computed as multiple of the
            full-cube RMS.
        scale_upper_limit : number, None
            Restrict the scales (in pixels) used by multi-scale clean to values
            less than this limit. If `None`, use all scales.
        """
        if scale_upper_limit is None:
            scales = self.scales
        else:
            scales = [s for s in self.scales if s < scale_upper_limit]
        assert len(scales) > 0
        dset, spw = self.dset, self.spw
        # Book-keeping before run
        imagename = self.nomask_imagebase
        log_post(':: Running clean ({0})'.format(imagename))
        # Calculate RMS values from off-line channels of dirty cube.
        threshold = format_rms(self.rms, sigma=sigma)
        # channel ranges: windowed or common coverage
        start, nchan = self.selected_start_nchan
        # run tclean
        delete_all_extensions(imagename)
        tclean(
            vis=dset.vis,
            imagename=imagename,
            field=dset.field,
            spw=self.spw_ids,
            specmode='cube',
            chanchunks=-1,
            outframe='lsrk',
            veltype='radio',
            restfreq=spw.restfreq,
            nchan=nchan,
            start=start,
            imsize=dset.imsize,
            cell=dset.cell,
            # gridder parameters
            gridder=dset.gridder,
            weighting=self.weighting_kind,
            robust=self.robust,
            perchanweightdensity=PERCHANWT,
            # deconvolver parameters
            deconvolver='multiscale',
            scales=scales,
            smallscalebias=self.smallscalebias,
            gain=self.gain,
            cyclefactor=self.cyclefactor,
            pblimit=dset.pblimit,
            niter=NITERMAX,
            threshold=threshold,
            interactive=False,
            parallel=self.parallel,
            # no masking parameters applied
            usemask='user',
        )
        self.cleanup(imagename)

    def make_seed_mask(self, sigma=5.0):
        imagename = self.nomask_imagebase
        if not os.path.exists(imagename+'.image'):
            log_post('-- File not found: {0}.image'.format(imagename))
            log_post('-- Has the unmasked call to tclean been run?')
            raise IOError
        make_multiscale_joint_mask(imagename, sigma=sigma,
                mask_ang_scales=self.mask_ang_scales)

    def clean_line(self, mask_method='seed+multithresh', sigma=2,
            ext=None):
        """
        Primary interface for calling `tclean` to deconvolve spectral windows.

        Parameters
        ----------
        mask_method : str
            Masking method to use in the deconvolution process. Available methods:
                'auto-multithresh' -- use auto-multithresh automated masking
                'seed+multithresh' -- generate initial mask from free clean
                    and then use auto-multithresh for automatic masking.
                'taper' -- use mask generated from separate tapering run
        sigma : number
            Threshold in standard deviations of the noise to clean down to within
            the clean-mask. An absolute RMS is calculated from off-line channels in
            a dirty cube. The same value is applied to all channels.
        ext : str
            String of the form '_EXT' appended to the end of the image name.
        """
        dset, spw = self.dset, self.spw
        # Book-keeping before run
        imagename = self.get_imagebase(ext=ext)
        log_post(':: Running clean ({0})'.format(imagename))
        # Calculate RMS values from off-line channels of dirty cube.
        threshold = format_rms(self.rms, sigma=sigma)
        # channel ranges: windowed or common coverage
        start, nchan = self.selected_start_nchan
        # masking method, auto-multithresh or mask from tapered run
        if mask_method == 'auto-multithresh' or mask_method == 'seed+multithresh':
            if dset.kind in ('joint', '12m'):
                mask_kwargs = self.autom_kwargs
            elif dset.kind == '7m':
                raise NotImplementedError
            else:
                raise ValueError('No parameters available for dset.kind: "{0}"'.format(dset.kind))
        elif mask_method == 'taper':
            # user supplied mask from seperate tapering run
            mask_kwargs = {
                    'usemask': 'user',
                    'mask': '{0}_taper.mask'.format(imagename),
                    }
        else:
            raise ValueError('Invalid mask_method: "{0}"'.format(mask_method))
        # run tclean
        delete_all_extensions(imagename)
        if mask_method == 'seed+multithresh':
            # copy summask to be used in-place with default extension
            # FIXME standard image doesn't appear to work with parallel=True
            mask_filen = self.nomask_imagebase + '.summask'
            shutil.copytree(mask_filen, imagename+'.mask')
        tclean(
            vis=dset.vis,
            imagename=imagename,
            field=dset.field,
            spw=self.spw_ids,
            specmode='cube',
            chanchunks=-1,
            outframe='lsrk',
            veltype='radio',
            restfreq=spw.restfreq,
            nchan=nchan,
            start=start,
            imsize=dset.imsize,
            cell=dset.cell,
            # gridder parameters
            gridder=dset.gridder,
            weighting=self.weighting_kind,
            robust=self.robust,
            perchanweightdensity=PERCHANWT,
            # deconvolver parameters
            deconvolver='multiscale',
            scales=self.scales,
            smallscalebias=self.smallscalebias,
            gain=self.gain,
            cyclefactor=self.cyclefactor,
            pblimit=dset.pblimit,
            niter=NITERMAX,
            threshold=threshold,
            interactive=False,
            parallel=self.parallel,
            # mask parameters
            **mask_kwargs
        )
        self.cleanup(imagename)

    def cleanup(self, imagename):
        delete_workdir(imagename)
        if self.parallel:
            concat_parallel_all_extensions(imagename)

    def postprocess(self, ext=None, make_fits=True):
        imagebase = self.get_imagebase(ext=ext)
        check_max_residual(imagebase, sigma=5.0)
        primary_beam_correct(imagebase)
        # smooth to common beam and export to FITS
        im_name = imagebase + '.image.pbcor'
        pb_name = imagebase + '.pb'
        cm_name = im_name   + '.common'
        smooth_cube_to_common_beam(im_name)
        copy_pb_mask(cm_name, pb_name)
        if make_fits:
            export_fits(cm_name)

    def run_pipeline(self, ext='clean'):
        self.make_dirty_cube()
        self.clean_line_nomask(sigma=4.5)
        self.make_seed_mask(sigma=5.0)
        self.clean_line(mask_method='seed+multithresh', ext=ext)
        self.postprocess(ext=ext, make_fits=True)


def make_all_line_dirty_cubes(dset, weighting=0.5, fullcube=True):
    setup = dset.setup
    log_post(':: Creating all dirty cubes for Setup-{0}'.format(setup))
    for spw in SPWS_BY_SETUP[setup].values():
        config = ImageConfig(dset, spw, weighting=weighting, fullcube=fullcube)
        config.make_dirty_cube()


def postproc_all_cleaned_images(dset):
    all_filen = glob('images/{0}/*.image'.format(dset.field))
    log_post(':: Post-processing image cubes to common beam')
    for imagename in all_filen:
        postproc_image_to_common(imagename)


def run_pipeline(field, setup=None, weightings=None, fullcube=True,
        do_cont=True):
    """
    Run all pipeline tasks over all setups, spectral windows, and
    uv-weightings.

    Parameters
    ----------
    field : str
    setup : (int, None)
        Setup number, if `None` run on all setups in sequence.
    weightings : Iterable, default (0.5,)
        List of uv-weightings to use in `tclean`, which may include the string
        "natural" or a number for the briggs robust parameter.
    fullcube : bool, default True
        Image the full spectral window or a small window around the SPWs
        defined rest frequency.
    do_cont : bool, default True
        Whether to image the continuum SPWs or not. 
    """
    weightings = WEIGHTINGS if weightings is None else weightings
    run_setups = (1, 2, 3) if setup is None else [setup]
    for i_setup in run_setups:
        dset = DataSet(field, setup=i_setup, kind='joint')
        for spw in SPWS_BY_SETUP[i_setup].values():
            if not do_cont:
                if spw.mol_name == 'cont':
                    continue
            log_post(':: Imaging Setup-{0} {1}'.format(setup, spw.label))
            for weighting in weightings:
                config = ImageConfig(dset, spw, fullcube=fullcube, weighting=weighting)
                config.run_pipeline()
    # TODO:
    #   - create moment maps
    #   - create diagnostic plots (channel maps and moments)


###############################################################################
# Tests and diagnostics
###############################################################################

def test_perchanwt():
    # mutate globals used within the scope of several functions
    global PERCHANWT
    global DIRTY_EXT
    for weighting in ('natural', 0.5):
        config = ImageConfig.from_name('CB68', '244.936GHz_CS', fullcube=True,
                weighting=weighting)
        for perchanwt in (True, False):
            PERCHANWT = perchanwt
            DIRTY_EXT = 'dirtyPCW{0}'.format(perchanwt)
            config.make_dirty_cube()
            config.clean_line_target(ext='PCW{0}'.format(perchanwt))
    # reset back to global defaults
    PERCHANWT = True
    DIRTY_EXT = 'dirty'


def test_rename_oldfiles(field, label=None, kind='joint', weighting=0.5):
    """
    Replace the rest frequency in each final image cube with the value
    given by Satoshi and then re-export the FITS cube.

    Parameters
    ----------
    field : str
        Field name, e.g., "CB68"
    label : (str, None)
        SPW label, e.g., "216.113GHz_DCOp". If `None` then apply to all
        Setups and SPWs for target field.
    kind : str
    weighting : (str, number)
    """
    log_post(':: Renaming final products from old frequency convention')
    if label is not None:
        spw_list = [ALL_SPWS[label]]
    else:
        spw_list = ALL_SPWS.values()
    for spw in spw_list:
        # setup new config instance and calculate frequencies
        dset = DataSet(field, setup=spw.setup, kind=kind)
        config = ImageConfig(dset, spw, weighting=weighting)
        new_base = config.get_imagebase(ext='clean')
        old_freq = '{0:.3f}GHz'.format(qa.convert(spw.ms_restfreq, 'GHz')['value'])
        old_base = 'images/{0}/{0}_{1}_{2}_{3}_{4}_clean'.format(
                config.dset.field, old_freq, spw.mol_name, kind, weighting,
        )
        old_name = '{0}.image.pbcor.common'.format(old_base)
        new_name = '{0}.image.pbcor.common'.format(new_base)
        if not os.path.exists(old_name):
            log_post('-- File not found, continuing: {0}'.format(old_name))
            continue
        # edit header rest frequency value
        restfreq_quantity = qa.quantity(spw.restfreq)  # in GHz
        restfreq_quantity = qa.convert(restfreq_quantity, 'Hz')
        log_post('-- Converting frequency for: {0}'.format(old_name))
        log_post('-- nu_old: {0:.12e}, nu_new: {1:.12e}'.format(
            float(spw.ms_restfreq.rstrip('Hz')),
            restfreq_quantity['value'],
        ))
        success = imhead(old_name, mode='put', hdkey='restfreq',
                hdvalue=restfreq_quantity)
        if not success:
            raise RuntimeError('-- Failed to modify header of: {0}'.format(old_name))
        # re-export FITS file
        export_fits(old_name)
        shutil.move(old_name+'.fits', new_name+'.fits')


###############################################################################
# Quality assurance plots
###############################################################################

def savefig(filen, dpi=300):
    for ext in ('png', 'pdf'):
        outfilen = os.path.join(PLOT_DIR, '{0}.{1}'.format(filen, ext))
        plt.savefig(outfilen, dpi=dpi)
    log_post('-- figure saved for: {0}'.format(filen))
    plt.close('all')


class CubeSet(object):
    def __init__(self, path, sigma=6):
        self.path = path
        self.sigma = sigma
        # Read in images and calculate noise
        stem = os.path.splitext(path)[0]
        self.stem = stem
        log_post('-- Reading image')
        self.image = self.get_chunk(stem+'.image')
        log_post('-- Reading residual')
        self.residual = self.get_chunk(stem+'.residual')
        log_post('-- Reading mask')
        self.mask = self.get_chunk(stem+'.mask')
        log_post('-- Reading primary beam')
        self.pb = self.get_chunk(stem+'.pb')
        self.shape = self.image.shape
        self.pix_width = self.shape[1]  # square images
        # calculate coordinate offset positions for ticks
        self.header = imhead(path, mode='list')
        self.pix_scale = abs(np.rad2deg(self.header['cdelt1']) * 60**2)  # arcsec
        # mask images on pb profile
        pb_mask = self.pb[0] < 0.2
        self.image[:,pb_mask] = np.nan
        self.residual[:,pb_mask] = np.nan
        # calculate noise
        self.rms = np.nanstd(self.image[:30])
        self.nchan = self.image.shape[0]
        # Identify channels with significant emisssion
        cut = self.image > sigma * self.rms
        self.good_mask = np.nansum(cut, axis=(1, 2)) > 0
        self.good_chan = np.argwhere(self.good_mask).flatten()
        self.ngood = self.good_chan.size

    @staticmethod
    def get_chunk(imagename):
        ia.open(imagename)
        # drop (degenerate) stokes axis
        chunk = ia.getchunk(dropdeg=True)
        ia.close()
        return chunk.transpose()

    def iter_planes(self):
        for ix in self.good_chan:
            im = self.image[ix]
            rs = self.residual[ix]
            ma = self.mask[ix]
            pb = self.pb[ix]
            yield im, rs, ma, pb, ix

    def calc_tick_loc(self, ang_tick=5):
        """
        Parameters
        ----------
        ang_tick : number
            Tick interval in arcsec.
        """
        ang_width = self.pix_scale * self.pix_width
        n_ticks = int(ang_width // ang_tick)
        n_ticks += 1 if n_ticks % 2 != 0 else 0
        pix_tick = ang_tick / self.pix_scale
        return np.arange(-n_ticks/2., n_ticks/2.+1) * pix_tick + self.pix_width / 2


def make_qa_plot(cset, kind='image', outfilen='qa_plot'):
    """
    Parameters
    ----------
    cset : CubeSet
    kind : str
        Kind of image data to plot:
            'image' : restored/cleaned image
            'residual' : residual image
    outfilen : str
    """
    log_post(':: Making QA plots for: {0}'.format(cset.stem))
    # configuration options specific to image or residual plots
    if kind == 'image':
        cmap = plt.cm.afmhot if kind == 'image' else plt.cm.RdBu
        mask_cont_color = 'cyan'
        vmin, vmax = None, None
    elif kind == 'residual':
        cmap = plt.cm.RdBu
        mask_cont_color = 'black'
        vmin, vmax = -5 * cset.rms, 5 * cset.rms
    else:
        raise ValueError('Invalid plot kind: "{0}"'.format(kind))
    # plot configuration options
    ncols = 5
    subplot_size = 1.3  # inch, good for five across on a 8.5x11 page
    # set NaN color for colormap, '0.2' for darker gray
    cmap.set_bad('0.5', 1.0)
    # determine the number of channels that need to be plotted
    nrows = cset.ngood // ncols + 1
    nempty = ncols - cset.ngood % ncols
    tick_pos = cset.calc_tick_loc(ang_tick=5)
    figsize = (ncols * subplot_size, nrows * subplot_size)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True,
            sharey=True, figsize=figsize)
    for ax, planes in zip(axes.flat, cset.iter_planes()):
        image, residual, mask, pbeam, chan_ix = planes
        data = image if kind == 'image' else residual
        # show colorscale of image data
        ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap, origin='lower')
        # show HPBW contour for primary beam data
        ax.contour(pbeam, levels=[0.5], colors=['0.5'],
                linestyles=['dashed'], linewidths=[0.4])
        # show contour for the clean mask
        if mask.any() > 0:
            ax.contour(mask, level=[0.5], colors=[mask_cont_color],
                    linestyles=['solid'], linewidths=[0.2])
        # show nice channel label in the top-right corner
        text = ax.annotate(str(chan_ix), (0.83, 0.91),
                xycoords='axes fraction', fontsize='xx-small')
        text.set_path_effects([
                path_effects.withStroke(linewidth=2, foreground='white')])
        # FIXME
        # show axis ticks get relative offsets in arcsec
        ax.set_xticks(tick_pos)
        ax.set_yticks(tick_pos)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    # zoom the window to crop out some of the PB NaNs around the edges
    ax.set_xlim(0.12*cset.pix_width, 0.88*cset.pix_width)
    ax.set_ylim(0.12*cset.pix_width, 0.88*cset.pix_width)
    plt.tight_layout(pad=0.8)
    # Hide unused plots
    for ax in axes.flat[-nempty:]:
        ax.set_visible(False)
    outfilen_ext = '{0}_{1}'.format(outfilen, kind)
    savefig(outfilen_ext)


def make_all_qa_plots(field, ext='clean'):
    image_paths = glob('{0}{1}/{1}_*_{2}.image'.format(IMAG_DIR, field, ext))
    for path in image_paths:
        cset = CubeSet(path)
        outfilen = '{0}_qa_plot'.format(cset.stem)
        make_qa_plot(cset, kind='image', outfilen=outfilen)
        make_qa_plot(cset, kind='residual', outfilen=outfilen)


###############################################################################
# User defined functions
###############################################################################

if __name__ == '__main__':
    # NOTE statements placed in this block will be executed on `execfile` in
    #      addition to when qsub is run.
    pass


