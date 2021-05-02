#!/usr/bin/env python
"""
=======================
FAUST Archival Pipeline
=======================
This module contains the implementation for the CASA spectral line imaging
pipeline for the archival products of the ALMA Large Program FAUST. Please
see the README file for further information and the documentation, which
may be found under the `docs/` directory or online at:

    https://faust-imaging.readthedocs.io

This script is to be run under CASA v5.6 and Python v2.7. To use the script
interactively, execute the module by running at the CASA IPython prompt:
``execfile('<PATH>/<TO>/faust_imaging.py')``.

Authors
-------
Brian Svoboda and Claire Chandler.

Copyright 2020, 2021 Brian Svoboda under the MIT License.
"""
from __future__ import (print_function, division)

import gc
import os
import sys
import shutil
import datetime
import warnings
from glob import glob
from copy import deepcopy
from collections import (OrderedDict, Iterable)

import numpy as np
from scipy import special
from matplotlib import pyplot as plt
from matplotlib import patheffects as path_effects

# This script is meant to be run under Python v2.7 in CASA v5.6. However,
# the documentation is generated using Sphinx and Python v3. In order for this
# script to execute its module level scope without crashing, a few global tasks
# in CASA must be mocked. Full CASA v5/v6 compatibility could be implemented if
# the user has CASA v6 installed with modules `casatasks` and `casatools`.
if sys.version_info < (3, 0, 0):
    # The script is run under CASA 5.6 Python 2.
    from cleanhelper import cleanhelper
else:
    # The script is run under a users Python v3 installation. This is only used
    # for generating the documentation.
    class Mock:
        is_mpi_enabled = True
        def convert(self, x, y):
            return {'unit': 'null', 'value': 1}
    qa = Mock()
    MPIEnvironment = Mock


# matplotlib configuration settings
plt.rc('font', size=10, family='serif')
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
plt.ioff()  # turn off interactive GUI


###############################################################################
# Global configuration options
###############################################################################

# Global paths. The `DATA_DIR` and `PROD_DIR` variables below are intended to
# be modified by the user for compatibility with their specific host system.
DATA_DIR = '/lustre/aoc/users/cchandle/FAUST/2018.1.01205.L/completed_SBs/'
PROD_DIR = '/lustre/aoc/users/bsvoboda/faust/faust_alma/data/'
IMAG_DIR = os.path.join(PROD_DIR, 'images/')
MOMA_DIR = os.path.join(PROD_DIR, 'moments/')
PLOT_DIR = os.path.join(PROD_DIR, 'plots/')

# Median-absolute-deviation conversion factor to Gaussian RMS
MAD_TO_RMS = 1 / (np.sqrt(2) * special.erfinv(0.5))  # approx. 1.4826
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
NOMASK_EXT = 'nomask'
# Default extension name for small image stamps used to calculate the ".sumwt"
# files.
TINYIMG_EXT = 'tinyimg'
# perchanelweightdensity parameter in tclean, changing to False slightly
# improves the resolution at the expense of sensitivity and uniform RMS.
PERCHANWT = False
# See guide for discussion for of automasking parameters:
#   https://casaguides.nrao.edu/index.php/Automasking_Guide
# These parameters have been tested to provide results that are
# frequently acceptable and are an overall compromise in masking compact
# and extended emission.
AUTOM_KWARGS = {
        'usemask': 'auto-multithresh',
        'noisethreshold': 5.0,
        'sidelobethreshold': 2.0,
        'lownoisethreshold': 1.5,
        'minbeamfrac': 0.001,
        'negativethreshold': 7.0,
        'growiterations': 1000,
        'dogrowprune': False,
        'fastnoise': False,
        'verbose': True,
}
# Default 7m values taken from the auto-masking CASA Guide above. Un-tested.
AUTOM_7M_KWARGS = {
        'usemask': 'auto-multithresh',
        'noisethreshold': 5.0,
        'sidelobethreshold': 1.25,
        'lownoisethreshold': 2.0,
        'minbeamfrac': 0.1,
        'negativethreshold': 10.0,
        'growiterations': 1000,
        'dogrowprune': False,
        'fastnoise': False,
        'verbose': True,
}


###############################################################################
# FAUST Target and SPW configuration
###############################################################################

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
    _ms_fmt = '{0}-Setup{1}/uid___*_target_lines_self_calibrated_continuum_subtracted_aligned.ms'
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
            values include: ('joint', '12m', '7m'). Note that only 12m data is
            availabe for Setup 3.
        """
        assert field in ALL_FIELD_NAMES
        assert setup in (1, 2, 3)
        if setup == 3 and kind != '12m':
            log_post('-- Reverting `kind` to "12m" for Setup 3', priority='WARN')
            kind = '12m'
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
    def ms_fmt(self):
        return os.path.join(DATA_DIR, self._ms_fmt)

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
        if self.setup == 3:
            return 'standard'
        elif self.setup in (1, 2) and self.kind == 'joint':
            return 'mosaic'
        elif self.setup in (1, 2) and self.kind in ('7m', '12m'):
            return 'standard'
        else:
            raise ValueError('Invalid setup, kind: ({0}, {1})'.format(self.setup, self.kind))

    @property
    def pblimit(self):
        return PBLIMIT if self.kind == 'joint' else -0.01

    @property
    def imsize(self):
        """
        Get the image size in pixels dynamically according to the lowest
        frequency SPW of the Setup and antenna diameter (i.e., 12m or 7m).
        The field size size is calculated as 10% larger than the full-width
        at the 10%-maximum point of the primary beam (approximated as
        $1.13 \lambda / D$).
        """
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

    def __init__(self, setup, restfreq, mol_name, name, ms_restfreq, spw_id,
            ot_name, nchan, chan_width, tot_bw):
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
            Channel width in TOPO frame.
            **units**: kHz.
        tot_bw : number
            Total bandwidth in TOPO frame.
            **units**: MHz
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

    def copy(self):
        return deepcopy(self)

    def with_chunk(self, chunk):
        """
        Parameters
        ----------
        chunk : ChunkConfig, None
        """
        if chunk is None:
            return self
        else:
            spw = self.copy()
            spw.mol_name = '{0}_chunk{1}'.format(self.mol_name, chunk.index)
            return spw


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

def log_post(msg, priority='INFO'):
    """
    Post a message to the CASA logger, logfile, and stdout/console.

    Parameters
    ----------
    msg : str
        Message to post.
    priority : str, default 'INFO'
        Priority level. Includes 'INFO', 'WARN', 'SEVERE'.
    """
    print(msg)
    casalog.post(msg, priority, 'faust_pipe')


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


def safely_remove_file(filen):
    if not os.path.exists(filen):
        log_post('-- File not found: {0}'.format(filen), priority='WARN')
        return
    log_post(':: Removing: {0}'.format(filen))
    if os.path.isdir(filen):
        # The target is a directory, i.e., a CASA image. Use `rmtables` to
        # safely release file locks within CASA on the internal tables.
        rmtables(filen)
        # The file will still exist if rmtables failed, such as for parallel
        # images that have a directory structure that `rmtables` does not
        # recognize.
        try:
            shutil.rmtree(filen)
            log_post('-- Hard delete: {0}'.format(filen))
        except OSError:
            pass
    else:
        # Target is a file and not a directory (e.g., FITS).
        os.remove(filen)


def if_exists_remove(filen):
    if os.path.exists(filen):
        safely_remove_file(filen)


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
        safely_remove_file(filen)


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
        msg = 'File not found: "{0}"'.format(new_filen)
        log_post(msg, priority='WARN')
        raise IOError(msg)
    rmtables(old_filen)
    # If a parallel image, then will not be recognized as a table in the call
    # to `rmtables`
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
    Create a contiguous image cube file from a 'parallel image' generated from
    a multiprocess `tclean` run with `parallel=True`. The function replaces the
    original parallel image with the contiguous image.
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
    ia.close()
    replace_existing_file_with_new(imagename, outfile)


def concat_parallel_all_extensions(imagebase):
    exts = ['image', 'mask', 'model', 'pb', 'psf', 'residual', 'sumwt', 'weight']
    for ext in exts:
        imagename = '{0}.{1}'.format(imagebase, ext)
        if os.path.exists(imagename):
            concat_parallel_image(imagename)


def get_freq_axis_from_image(imagename):
    """
    Compute the frequency axis from an image's header coordinate system
    information.

    Parameters
    ----------
    imagename : str

    Returns
    -------
    ndarray
        Frequency axis.
    str
        Frequency axis unit (likely Hz).
    """
    # Get coordinate system information.
    ia.open(imagename)
    csys = ia.coordsys()
    csys_dict = csys.torecord()
    shape = ia.shape()
    ia.close()
    ia.done()
    # Validate the axes and get the number of channels.
    try:
        spec_index = csys.findaxisbyname('Frequency')
    except RuntimeError:
        raise ValueError('Image has no frequency axis, is in velocity?')
    nchan = shape[spec_index]
    csys.done()
    # Compute the frequency axis indices to frequencies.
    spec_sys = csys_dict['spectral2']
    sunit = spec_sys['unit']
    crpix = spec_sys['wcs']['crpix']
    crval = spec_sys['wcs']['crval']
    cdelt = spec_sys['wcs']['cdelt']
    faxis = (np.arange(nchan) + crpix) * cdelt + crval
    return faxis, sunit


def calc_common_coverage_range(imagebase):
    """
    Calculate the common frequency coverage amongst a set of MSs/EBs. The
    `.sumwt` image extension is used from an existing "dirty" cube.

    Parameters
    ----------
    imagebase : str
        Image name without the ".sumwt" extension.

    Returns
    -------
    str
        The start frequency in the LSRK frame.
    int
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
    faxis, sunit = get_freq_axis_from_image(sumwt_filen)
    start_val = faxis[ix_lo:ix_hi].min()
    start = '{0}{1}'.format(start_val, sunit)
    log_post('-- start frequency: {0}'.format(start))
    log_post('-- nchan: {0}'.format(nchan))
    return start, nchan


def calc_chunk_freqs(imagebase, nchunks=1):
    """
    Divide the spectral axis of an image into continuous chunks in frequency.
    Returned frequencies specify the channel center frequency.

    Parameters
    ----------
    imagebase : str
        Image name base without the ".sumwt" extension.
    nchunks : int
        Number of contiguous chunks to divide the frequency range into.

    Returns
    -------
    [[str, int], ...]
        (str) LSRK frequency of the bin left-edge of the first channel.
        **units**: Hz.

        (int) Number of channels in the chunk.
    """
    nchunks = int(nchunks)
    assert nchunks >= 1
    freq_start, nchan_total = calc_common_coverage_range(imagebase)
    assert nchunks <= nchan_total
    # Convert start frequency unit string into a float.
    spec_unit = 'Hz'
    freq_start = qa.convert(freq_start, spec_unit)['value']
    # Get frequency channel width.
    imagename = '{0}.sumwt'.format(imagebase)
    ia.open(imagename)
    csys = ia.coordsys()
    csys_dict = csys.torecord()
    ia.close()
    ia.done()
    csys.done()
    # Ensure that the delta and start frequency have the same units.
    cdelt = csys_dict['spectral2']['wcs']['cdelt']
    cunit = csys_dict['spectral2']['unit']
    cdelt = qa.convert(qa.quantity(cdelt, cunit), spec_unit)['value']
    cdelt = abs(cdelt)
    # Calculate the number of channels per chunk. Distribute the channels the
    # channel numbers that don't evenly divide uniformly among the remaining.
    if nchan_total == nchunks:
        chan_per_chunk = np.ones(nchunks, dtype=int)
    else:
        chan_per_chunk = np.zeros(nchunks, dtype=int)
        chan_per_chunk[:] = nchan_total // nchunks
        remainder = nchan_total % nchunks
        chan_per_chunk[-remainder:] += 1
        assert chan_per_chunk.sum() == nchan_total
    # Calculate the starting frequencies by iteratively adding the chunk
    # bandwidth to the last starting frequency.
    chunks = []
    freq = freq_start
    for n in chan_per_chunk:
        freq_withunit = qa.tos(qa.quantity(freq, spec_unit))
        chunks.append([freq_withunit, n])
        freq += n * cdelt
    return chunks


def primary_beam_correct(imagebase):
    log_post(':: Primary beam correcting image: {0}'.format(imagebase))
    impbcor(
            imagename=imagebase+'.image',
            pbimage=imagebase+'.pb',
            outfile=imagebase+'.image.pbcor',
            overwrite=True,
    )


def hanning_smooth_image(imagename, overwrite=True):
    """
    Hanning smooth an image. Produces a new file with the extension ".hanning".

    Parameters
    ----------
    filen : str
    overwrite : bool
    """
    log_post(':: Creating Hanning smoothed image')
    outname = '{0}.hanning'.format(imagename)
    ia.open(imagename)
    ia.hanning(
            outname, drop=False, dmethod='mean', overwrite=overwrite,
    ).done()
    ia.done()
    ia.close()


def smooth_cube_to_common_beam(imagename, beam=None):
    """
    Use the ``imsmooth`` task to convert a cube with per-plane beams into one
    with a single common beam.

    Parameters
    ----------
    imagename : str
        Name of the image, including extension such as ".image".
    beam : dict
        Gaussian beam parameter dictionary (see ``imsmooth`` help).
    """
    log_post(':: Smoothing cube to common beam scale')
    outfile = imagename + '.common'
    if beam is None:
        imsmooth(imagename=imagename, kernel='commonbeam', outfile=outfile,
                overwrite=True)
    else:
        imsmooth(imagename=imagename, kernel='gaussian', outfile=outfile,
                beam=beam, targetres=True, overwrite=True)


def copy_pb_mask(imagename, pbimage):
    """
    Copy primary beam T/F mask back into the image after being removed by a
    task such as ``imsmooth``.
    """
    log_post(':: Copying primary beam mask into: {0}'.format(imagename))
    makemask(mode='copy', inpimage=pbimage, inpmask=pbimage+':mask0',
            output=imagename+':mask0', overwrite=True)


def calc_rms_from_image(imagename, chan_expr=None):
    """
    Calculate RMS values from the scaled MAD of all channels (unless given a
    range) for the given full imagename.

    Parameters
    ----------
    imagename : str
        CASA Image name, e.g., "244.936GHz_CS_joint_0.5_dirty.image"
    chan_expr : str, None
        Channel selection expression passed to `imstat` parameter `chans`.
        If None, the full channel range is used.

    Returns
    -------
    float
        RMS
    """
    if not os.path.exists(imagename):
        msg = '-- File or directory not found: "{0}"'.format(imagename)
        log_post(msg, priority='WARN')
        raise IOError(msg)
    # run imstat over the cube and extract the calculated RMS
    stats = imstat(imagename=imagename, axes=[0, 1], chans=chan_expr)
    mad_arr = stats['medabsdevmed']
    rms = MAD_TO_RMS * np.nanmedian(mad_arr)
    log_post('-- RMS: {0}'.format(rms))
    return rms


def create_mask_from_threshold(infile, outfile, thresh, overwrite=True):
    """
    Create a 1/0 mask according to whether pixel values in the input image are
    or are not above ``thresh`` in Jy.

    Parameters
    ----------
    infile : str
    outfile : str
    thresh : number
        **units**: Jy
    overwrite : bool, default True
    """
    if overwrite:
        if_exists_remove(outfile)
    immath(imagename=infile, mode='evalexpr', outfile=outfile,
            expr='iif(IM0>{0},1,0)'.format(thresh))


def effective_beamwidth_from_image(imagename):
    """
    Median of the geometric-mean beamwidths in an image with perplane beams.
    If the image has a single plane or single restoring beam, then the
    geometric mean of the major & minor axes is returned.

    Parameters
    ----------
    imagename : str
    """
    if not os.path.exists(imagename):
        msg = '-- File not found to compute beamwidth: "{0}"'.format(imagename)
        log_post(msg, priority='WARN')
        raise IOError(msg)
    header = imhead(imagename)
    if 'perplanebeams' in header:
        geom_mean_widths = np.sqrt([
                v['*0']['major']['value'] * v['*0']['minor']['value']
                for v in header['perplanebeams']['beams'].values()
        ])
        return np.median(geom_mean_widths)
    elif 'restoringbeam' in header:
        beam = header['restoringbeam']
        geom_mean_width = np.sqrt(beam['major']['value'] * beam['minor']['value'])
        return geom_mean_width
    else:
        raise ValueError('Invalid beam information in: {0}'.format(imagename))


def make_multiscale_joint_mask(imagename, rms, sigma=5.0, mask_ang_scales=(0, 1, 3)):
    """
    Create a joint mask from multiple smoothed copies of an image. A file
    with the `.summask` extension is created from the union of an RMS threshold
    applied to each image. The RMS is re-computed for each smoothed cube.

    Parameters
    ----------
    imagename : str
        Image name base without the extension.
    rms : number
        RMS value of the unsmoothed image.
        **units**: Jy
    sigma : number, Iterable(number)
        Multiple of the RMS to threshold each image. If passed as a number,
        the same threshold is applied to all images (whether unsmoothed or
        smoothed). If passed as a list-like, then the sigma value is used for
        the image with the corresponding scale in `mask_ang_scales`.
    mask_ang_scales : Iterable(number)
        Gaussian kernel FWHM in units of arcseconds to convolve each image
        with (note: not in pixel units).

    Returns
    -------
        Files are written for each smoothing scale:
            smoothed image: '<IMG>_smooth{.3f}.image' (excluding scale=0)
            masked image:   '<IMG>_smooth{.3f}.mask'
        and the joint or unioned mask file across scales is written to:
            joint mask:     '<IMG>.summask'
    """
    if isinstance(sigma, Iterable):
        assert len(sigma) == len(mask_ang_scales)
        thresholds = [s * rms for s in sigma]
    else:
        thresholds = len(mask_ang_scales) * [sigma * rms]
    log_post(':: Creating threshold based mask')
    n_scales = len(mask_ang_scales)
    # Set file names and paths
    original_image = '{0}.image'.format(imagename)
    pb_image = '{0}.pb'.format(imagename)
    base_names = [
            '{0}_smooth{1:.3f}'.format(imagename, s)
            for s in mask_ang_scales
    ]
    mask_files = [s+'.mask' for s in base_names]
    image_files = [s+'.image' for s in base_names]
    # Smooth and threshold the image for each angular scale
    beamwidth = effective_beamwidth_from_image(original_image)
    for thresh, scale, smooth_image, mask_image in zip(
            thresholds, mask_ang_scales, image_files, mask_files):
        assert scale >= 0
        if scale == 0:
            # For a smoothing scale of zero, use the original image.
            create_mask_from_threshold(original_image, mask_image, thresh)
            continue
        # Smooth the image with a Gaussian kernel of the given scale.
        if_exists_remove(smooth_image)
        imsmooth(imagename=original_image, kernel='gauss',
                major='{0}arcsec'.format(scale),
                minor='{0}arcsec'.format(scale), pa='0deg', targetres=True,
                outfile=smooth_image)
        # PB mask is erased after smoothing, so copy it back in from the PB
        copy_pb_mask(smooth_image, pb_image)
        # Threshold the image to create a mask. To account for the different
        # RMS in the smoothed maps (implicitly in Jy/beam), the threshold needs
        # to be scaled according to the relative (target) beam size.
        beam_factor = max(1, scale / beamwidth)
        create_mask_from_threshold(smooth_image, mask_image, beam_factor*thresh)
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

class ChunkConfig(object):
    def __init__(self, start, nchan, index):
        """
        Parameters
        ----------
        start : str
        nchan : int
        index : int
        """
        assert qa.isquantity(start)
        self.start = start
        self.nchan = nchan
        self.index = int(index)


class ImageConfig(object):
    _scales = [0, 15, 45, 135]  # pixels; point, 1.5, 4.5 beam hpbw's (10 pix)
    _mask_ang_scales = [0, 1, 3]  # arcsec
    _autom_kwargs = AUTOM_KWARGS
    _autom_7m_kwargs = AUTOM_7M_KWARGS
    _rms = None
    smallscalebias = -1.0
    gain = 0.05
    cyclefactor = 2.0
    parallel = MPIEnvironment().is_mpi_enabled

    def __init__(self, dset, spw, fullcube=True, weighting=0.5, chunk=None):
        """
        Configuration object for `tclean` related custom tasks. Parameters
        specify deconvolution and image-cube properties such as the weighting.
        See the docstring for :meth:`from_name` for more complete
        parameter descriptions.

        Parameters
        ----------
        dset : DataSet
        spw : Spw
        fullcube : bool
        weighting : str, number
        chunk : ChunkConfig, None

        Attributes
        ----------
        scales : Iterable
            List of scales in pixels to be used for multi-scale clean. Will be
            passed to the ``scales`` keyword argument in `tclean`.
        mask_ang_scales : Iterable
            FWHM sizes in arcsec for Gaussian smoothing kernels used in
            :meth:`faust_imaging.ImageConfig.clean_line_nomask`.
        autom_kwargs : dict
            Auto-multithresh keyword arguments passed to `tclean`.
        smallscalebias : number
        gain : number
        cyclefactor : number
        parallel : bool
            Is MPI enabled in CASA?
        """
        assert dset.setup == spw.setup
        self.dset = dset
        self.spw = spw.with_chunk(chunk)
        self.fullcube = fullcube
        self.chunk = chunk
        if not fullcube and chunk is not None:
            raise ValueError("fullcube must be `True` when imaging a chunk.")
        # weighting method, must be valid for tclean
        robust, weighting_kind = format_cube_robust(weighting)
        self.weighting = weighting
        self.robust = robust
        self.weighting_kind = weighting_kind
        # copy mutable class attributes to avoid global mutation in all instances
        self.scales = deepcopy(self._scales)
        self.mask_ang_scales = deepcopy(self._mask_ang_scales)
        # set auto-multithresh keyword arguments for different array configurations
        if dset.kind in ('joint', '12m'):
            self.autom_kwargs = deepcopy(self._autom_kwargs)
        elif dset.kind == '7m':
            self.autom_kwargs = deepcopy(self._autom_7m_kwargs)
        else:
            raise ValueError('No default auto-multithresh parameters available for dset.kind: "{0}"'.format(dset.kind))

    @classmethod
    def from_name(cls, field, label, kind='joint', **kwargs):
        """
        Create an `ImageConfig` object directly from the field name and SPW
        label without creating instances of auxillary class instances.

        Parameters
        ----------
        field : str
            Field name, e.g., "CB68"
        label : str, None
            SPW label, e.g., "216.113GHz_DCOp". If `None` then apply to all
            Setups and SPWs for target field.
        kind : str
            Datset descriptor for which measurement set files to use. Valid
            values include: ('joint', '12m', '7m').
        fullcube : bool
            Image the full spectral window or a small window around the SPWs
            defined rest frequency.
        weighting : str, number
            Either a string for a weighting method recognized by `tclean`,
            such as 'natural' or 'uniform'; or a number for a Briggs weighting
            robust value.
        chunk : ChunkConfig, None
            Chunk configuration instance for dividing the frequency axis into
            individual, smaller memory foot-print intervals in frequency. The
            default value of `None` will apply no chunking along the frequency
            axis.
        """
        spw = ALL_SPWS[label]
        dset = DataSet(field, setup=spw.setup, kind=kind)
        return cls(dset, spw, **kwargs)

    @property
    def is_chunked(self):
        return self.chunk is not None

    @property
    def dirty_imagebase(self):
        return self.get_imagebase(ext=DIRTY_EXT)

    @property
    def nomask_imagebase(self):
        return self.get_imagebase(ext=NOMASK_EXT)

    @property
    def tinyimg_imagebase(self):
        return self.get_imagebase(ext=TINYIMG_EXT)

    @property
    def rms(self):
        """
        Image RMS. By default, determined using `imstat` over the entire dirty
        cube. The RMS can also be manually set.
        """
        if self._rms is None:
            imagename = '{0}.image'.format(self.dirty_imagebase)
            if not os.path.exists(imagename):
                msg = '-- File not found to compute RMS: "{0}"'.format(imagename)
                log_post(msg, priority='WARN')
                log_post('-- Has the `.make_dirty_cube` been run?')
                raise IOError(msg)
            rms = calc_rms_from_image(imagename)
            return rms
        else:
            return self._rms

    @rms.setter
    def rms(self, rms):
        self._rms = rms

    @property
    def selected_start_nchan(self):
        if self.fullcube and not self.is_chunked:
            # Restrict coverage to common coverage cross EBs.
            return calc_common_coverage_range(self.dirty_imagebase)
        elif self.is_chunked:
            return self.chunk.start, self.chunk.nchan
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

    def duplicate_into_chunks(self, nchunks=None):
        """
        Duplicate the image configuration into multiple versions set to be chunked
        in frequency. This eases memory requirements because each smaller image
        cubes may be processed independently. A full-bandwidth ".sumwt" file must
        exist to compute the starting frequencies. If neither
        :meth:`faust_imaging.ImageConfig.make_dirty_cube` nor
        :meth:`faust_imaging.ImageConfig.make_tiny_image_stamp` have been executed,
        the latter will be run before calculating the chunk frequencies.

        Parameters
        ----------
        nchunks : int, None
            Number of chunks to create. If unset then the number of chunks is
            chosen using a heuristic.

        Returns
        -------
        ChunkedConfigSet
            Configuration object encapsulating `ImageConfig` instances with
            properties set for chunking.
        """
        if self.is_chunked:
            raise ValueError("Config is already chunked: {0}".format(self.chunk))
        if nchunks is None:
            is_s12 = self.spw.setup in (1, 2)
            is_fdm = self.spw.mol_name != 'cont'
            nchunk_map = {
                    #  S12?   FDM?
                    ( True,  True):   4,  # S1/S2 FDM
                    ( True, False):  40,  # S1/S2 TDM
                    (False,  True):  20,  # S3    FDM
                    (False, False):  80,  # S3    TDM
            }
            nchunks = nchunk_map[(is_s12, is_fdm)]
        assert isinstance(nchunks, int) and nchunks > 0
        # A full-bandwidth cube must exist to compute the frequency ranges. Use
        # the full dirty image if it exists, otherwise use or make a tiny image
        # stamp for those purposes.
        if os.path.exists(self.dirty_imagebase+'.sumwt'):
            imagebase = self.dirty_imagebase
        elif os.path.exists(self.tinyimg_imagebase+'.sumwt'):
            imagebase = self.tinyimg_imagebase
        else:
            self.make_tiny_image_stamp()
            imagebase = self.tinyimg_imagebase
        # Create chunk configuration instances from the file-path to the image
        # products including the ".sumwt" file.
        chunks = [
                ChunkConfig(start, nchan, i)
                for i, (start, nchan) in enumerate(calc_chunk_freqs(
                    imagebase, nchunks=nchunks))
        ]
        chunked_configs = []
        # Copy and mutate the existing instance into modified chunked versions.
        # NOTE If a more sophisticated configuration is ever added to the
        #      initialization, then this simple over-writing of the attributes
        #      may leave instances improperly initialized.
        for chunk in chunks:
            log_post('-- Chunk start frequency: ({0}, {1})'.format(chunk.index, chunk.start))
            config = deepcopy(self)
            config.spw = config.spw.with_chunk(chunk)
            config.chunk = chunk
            config.chunked = True
            chunked_configs.append(config)
        return ChunkedConfigSet(chunked_configs)

    def get_imagebase(self, ext=None):
        """
        Get the relative path name for an image basename following

            "images/<FIELD>/<LABEL>_<ARRAY>_<WEIGHTING>[_<EXT>]"

        Note that this does not include image extensions such as ".image"

        Parameters
        ----------
        ext : str, Iterable, object
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

    def cleanup(self, imagename):
        """
        Delete the ".workdirectory" folders that are occassionally generated
        and not removed in mpicasa. Also concatenate all image products into
        contiguous "serial" products using `ia.imageconcat`.
        """
        delete_workdir(imagename)
        if self.parallel:
            concat_parallel_all_extensions(imagename)

    def remove_all_files(self, confirm=True):
        """
        Attempt to safely remove all files matching the image basename
        returned by :meth:`faust_imaging.ImageConfig.get_imagebase`.

        Parameters
        ----------
        confirm : bool
            Require confirmation before removing files.
        """
        imagebase = self.get_imagebase()
        matched_files = glob('{0}*'.format(imagebase))
        matched_files.sort()
        if confirm:
            for filen in matched_files:
                print('  ', filen)
            response = raw_input('Delete the above files? [y/N] ')
            if response.lower() != 'y':
                return
        for filen in matched_files:
            safely_remove_file(filen)

    def make_tiny_image_stamp(self):
        """
        Image a small field in order to produce the ".sumwt" file without
        producing a full set of dirty image products (which may be very large).
        The default suffix set by ``TINYIMG_EXT`` is "_tinyimg".
        """
        dset, spw = self.dset, self.spw
        imagename = self.tinyimg_imagebase
        log_post(':: Creating cube for sumwt ({0})'.format(imagename))
        # Create a dirty cube using the full-bandpass.
        niter = 0
        start = None
        nchan = -1
        # For this purpose the cell size does not need to be as aggressively
        # oversampled. Field sizes of about 10 HPBW appear to be the
        # approximate mininum required to create the PSF file and compute the
        # weights. Use a bit larger to be safe.
        oversampling_factor = 3
        cell = qa.tos(qa.quantity(dset.target.res/oversampling_factor, 'arcsec'))
        imwidth = cleanhelper.getOptimumSize(int(15*oversampling_factor))
        imsize = [imwidth, imwidth]
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
            imsize=imsize,
            cell=cell,
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
        # Validate that the sumwt file contains useful data. If tclean fails to
        # create the PSF then the sumwt file will contain all zeros.
        sumwt_filen = '{0}.sumwt'.format(imagename)
        ia.open(sumwt_filen)
        data = ia.getchunk()
        ia.close()
        ia.done()
        if np.all(data == 0):
            log_post('-- Failed to make PSF file, sumwt is all zeros: "{0}"'.format(sumwt_filen))
            raise RuntimeError

    def make_dirty_cube(self):
        """
        Generate a dirty image cube and associated `tclean` products.
        Run with equivalent parameters as `clean_line` but with `niter=0`.
        The default suffix set by ``DIRTY_EXT`` is "_dirty".
        """
        dset, spw = self.dset, self.spw
        # Book-keeping before run
        imagename = self.dirty_imagebase
        log_post(':: Creating dirty cube ({0})'.format(imagename))
        # To create dirty cube simply perform zero iterations.
        niter = 0
        # When imaging the full cube without chunking, image the full spectral
        # coverage over all array configurations (nchan=-1) but use the
        # pre-computed range when imaging a chunk.
        if self.is_chunked:
            start, nchan = self.selected_start_nchan
        else:
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
        Deconvolve the image without using a mask to a depth of `sigma` and
        excluding multiscale size scales beyond `scale_upper_limit`.

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
        """
        Create a mask based on a significance cut applied to the restored image
        of the unmasked run generated from `clean_line_nomask`. This mask is
        used to "seed" the mask generated using auto-multithresh in `clean_line`
        with `mask_method='seed+multithresh'`.

        Parameters
        ----------
        sigma : number
            Limit to threhsold the restored image on.
        """
        imagename = self.nomask_imagebase
        if not os.path.exists(imagename+'.image'):
            msg = '-- File not found: "{0}.image"'.format(imagename)
            log_post(msg, priority='WARN')
            log_post('-- Has the unmasked call to tclean been run?')
            raise IOError(msg)
        make_multiscale_joint_mask(imagename, self.rms, sigma=sigma,
                mask_ang_scales=self.mask_ang_scales)

    def clean_line(self, mask_method='seed+multithresh', sigma=3.0,
            ext=None, restart=False, interactive=False):
        """
        Primary interface for calling `tclean` to deconvolve spectral windows.
        Multiple masking methods exist to automatically generate the clean
        masks used. The preferred method is `mask_method='seed+multithresh'`,
        which requires that `.clean_line_nomask` is run first.

        Parameters
        ----------
        mask_method : str
            Masking method to use in the deconvolution process. Available methods

                ``"auto-multithresh"`` use auto-multithresh automated masking.

                ``"seed+multithresh"`` generate initial mask from free clean
                and then use auto-multithresh for automatic masking.

                ``"fixed"`` use the existing clean mask but keep it fixed and
                do not apply auto-multithresh.

                ``"taper"`` use mask generated from a separate tapered run.
                Generating the tapared masks requires re-implementation.

        sigma : number
            Threshold in standard deviations of the noise to clean down to within
            the clean-mask. An absolute RMS is calculated from off-line channels in
            a dirty cube. The same value is applied to all channels.
        ext : str
            String of the form '_EXT' appended to the end of the image name.
        restart : bool, default False
            Restart using the existing model and mask files.
        interactive : bool, default False
            Begin tclean in interactive mode with `interactive=True`. This may
            be useful for touching-up some channels with particularly difficult
            to clean extended emission.
        """
        dset, spw = self.dset, self.spw
        # Book-keeping before run
        imagename = self.get_imagebase(ext=ext)
        log_post(':: Running clean ({0})'.format(imagename))
        # Calculate RMS values from off-line channels of dirty cube.
        threshold = format_rms(self.rms, sigma=sigma)
        # channel ranges: windowed or common coverage
        start, nchan = self.selected_start_nchan
        # masking method, auto-multithresh or user masks
        if mask_method == 'auto-multithresh' or mask_method == 'seed+multithresh':
            mask_kwargs = self.autom_kwargs
        elif mask_method == 'existing':
            default_mask = '{0}.mask'.format(imagename)
            if not os.path.exists(default_mask):
                raise RuntimeError('Mask does not exist for mask_method="existing": {0}'.format(default_mask))
            mask_kwargs = {
                    'usemask': 'user',
                    'mask': default_mask,
                    }
        elif mask_method == 'taper':
            # user supplied mask from seperate tapering run
            mask_kwargs = {
                    'usemask': 'user',
                    'mask': '{0}_taper.mask'.format(imagename),
                    }
        else:
            raise ValueError('Invalid mask_method: "{0}"'.format(mask_method))
        # run tclean
        if restart:
            safely_remove_file('{0}.image'.format(imagename))
        else:
            delete_all_extensions(imagename)
        if not restart and mask_method == 'seed+multithresh':
            if self.parallel:
                # FIXME The mask generated is a "serial" image which is not
                # read correctly, and subsequently ignored, by tclean run in
                # parallel. It appears the image must be converted to a
                # "parallel image" first somehow.
                raise NotImplementedError('mask_method="seed+multithresh" unsupported with parallel=True')
            # copy summask to be used in-place with default extension
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
            # general run properties
            interactive=interactive,
            parallel=self.parallel,
            # mask parameters
            **mask_kwargs
        )
        self.cleanup(imagename)

    def clean_line_interactive_restart(self, **kwargs):
        self.clean_line(restart=True, interactive=True, **kwargs)

    def postprocess(self, ext=None, make_fits=True, make_hanning=True,
            beam=None):
        """
        Post-process image products. The maximum residual is logged to the
        CASA log file and STDOUT, the image is primary beam corrected, smoothed
        to a common beam, and the final FITS file is exported.

        Parameters
        ----------
        ext : str
        make_fits : bool
        make_hanning : bool
        beam : dict, None
            Gaussian beam parameter dictionary to use for common beam smoothing.
            If left as None, then the common beam of the image is used.
        """
        imagebase = self.get_imagebase(ext=ext)
        check_max_residual(imagebase, sigma=5.5)
        primary_beam_correct(imagebase)
        # smooth to common beam and export to FITS
        pb_name = imagebase + '.pb'
        for im_ext in ('.image', '.image.pbcor'):
            im_name = imagebase + im_ext
            cm_name = im_name + '.common'
            smooth_cube_to_common_beam(im_name, beam=beam)
            copy_pb_mask(cm_name, pb_name)
        if make_hanning:
            hanning_smooth_image(imagebase+'.image.common')
        if make_fits:
            image_to_export = imagebase + '.image.pbcor.common'
            export_fits(image_to_export)

    def run_pipeline_tasks(self, ext='clean', nomask_sigma=4.5,
            seedmask_sigma=5.0, clean_sigma=3.0):
        """
        Helper method to run the pipeline tasks in sequence with reasonable
        default parameters. See :meth:`run_pipeline` for further description.

        Parameters
        ----------
        ext : str
        nomask_sigma : number
        seedmask_sigma : number, Iterable(number)
        clean_sigma : number
        """
        self.make_dirty_cube()
        self.clean_line_nomask(sigma=nomask_sigma)
        self.make_seed_mask(sigma=seedmask_sigma)
        self.clean_line(mask_method='seed+multithresh', ext=ext, sigma=clean_sigma)
        # Unnecessary to post-process a chunked image.
        if not self.is_chunked:
            self.postprocess(ext=ext, make_fits=True)

    def run_pipeline(self, ext='clean', use_chunking=True, nchunks=None,
            nomask_sigma=4.5, seedmask_sigma=5.0, clean_sigma=3.0):
        """
        Run all pipeline tasks with default parameters. Custom recipes should
        call the individual methods in sequence. The final pipeline product
        will be named of the form:

            "<FILENAME>_<EXT>.image.pbcor.common.fits"

        or if ``ext=None``:

            "<FILENAME>.image.pbcor.common.fits"

        For convenience, some significance thresholds ("sigma") may also be set
        for this method. For more customized pipeline procedures, calling the
        individual pipeline methods is recommended.

        Parameters
        ----------
        ext : str
            Extension name for final, deconvolved image products.
        use_chunking : bool
            If True then chunk the image into separate frequency intervals,
            process individually in serial, and then concatenate the final
            results.
        nchunks : int, None
            Number of approximately uniform frequency intervals to chunk image
            products into. If `None` then the number of chunks is chosen
            by heuristic.
        nomask_sigma : number
            Global clean threshold for the un-masked clean run.
        seedmask_sigma : number, Iterable(number)
            Significance with which to threshold the un-masked clean run to
            create the seed mask. If a list-like is passed, the corresponding
            significance is used for each smoothing scale (scales are set by
            :attr:`.mask_ang_scales`).
        clean_sigma : number
            Global clean threshold of the final clean run.
        """
        sigma_kwargs = {'nomask_sigma': nomask_sigma, 'seedmask_sigma':
                seedmask_sigma, 'clean_sigma': clean_sigma}
        # If already chunked or explicitly selected to not use chunking, then
        # do not perform additional chunking and simply run the pipeline tasks.
        if self.is_chunked or not use_chunking:
            self.run_pipeline_tasks(ext=ext, **sigma_kwargs)
        else:
            chunked_configs = self.duplicate_into_chunks(nchunks=nchunks)
            for config in chunked_configs:
                config.run_pipeline_tasks(ext=ext, **sigma_kwargs)
            chunked_configs.postprocess(ext=ext)


class ChunkedConfigSet(object):
    def __init__(self, configs):
        """
        Parameters
        ----------
        configs : Iterable(ImageConfig)
            List of chunked ImageConfig instances, such as generated by
            :meth:`ImageConfig.duplicate_into_chunks`.

        Attributes
        ----------
        configs : Iterable(ImageConfig)
        nchunks : int
        nchan_total : int
            Total number of channels summing across all chunks.
        """
        assert isinstance(configs, Iterable)
        self.configs = configs
        self.nchunks = len(configs)
        self.nchan_total = sum(c.chunk.nchan for c in configs)

    def __iter__(self):
        for config in self.configs:
            yield config

    def __getitem__(self, key):
        return self.configs[key]

    def get_imagebase(self, ext=None):
        # Remove "_chunk0" in the first chunk to get the outfile base name.
        first_config = self.configs[0]
        assert '_chunk0' in first_config.dirty_imagebase
        imagebase = first_config.get_imagebase(ext=ext).replace('_chunk0', '')
        return imagebase

    def get_chunk_from_channel(self, ix):
        """
        Parameters
        ----------
        ix : int
            Channel index number in the full, concatenated cube.

        Returns
        -------
        ImageConfig
        """
        assert isinstance(ix, int)
        last_valid_index = self.nchan_total - 1
        assert 0 <= ix <= last_valid_index
        start = 0
        for config in self:
            nchan = config.chunk.nchan
            if start <= ix < start+nchan:
                return config
            start += nchan
        raise RuntimeError('Could not find channel index in chunks: {0}'.format(ix))

    def get_common_beam(self, ext='clean'):
        """
        Calculate the common beam (largest area) among all chunked images.

        Parameters
        ----------
        ext : str
            Extension name suffix, e.g. 'clean', 'nomask', etc.

        Returns
        -------
        (str, str, str)
            Beam parameters formatted as strings with unit labels: major axis
            (arcsec), minor axis (arcsec), position angle (deg).

        Notes
        -----
        A fudge factor of 0.3% is added to the beam size to account for
        discrepencies in how the common beam is calculated by CASA over the
        full cube compared to the largest-area common beam among the chunks.
        For typical ~0.35as HPBWs in the FAUST program, this multiplicatie
        factor corresponds to ~1 mas.
        """
        # Iteratively identify the largest beam among the common beams of all
        # the chunked images.
        max_area = 0
        common = None
        for config in self.configs:
            imagebase = config.get_imagebase(ext=ext)
            imagename = '{0}.image'.format(imagebase)
            ia.open(imagename)
            beam = ia.commonbeam()
            ia.close()
            ia.done()
            area = beam['major']['value'] * beam['minor']['value']
            assert not np.isnan(area)
            if area > max_area:
                max_area = area
                common = beam
        assert common is not None
        # Apply fudge factor
        epsilon = 1.003
        common['major']['value'] *= epsilon
        common['minor']['value'] *= epsilon
        return common

    def concat_cubes(self, ext='clean', im_exts=None):
        """
        Use ``ia.imageconcat`` to concatenate chunked image files into single
        a cube with contiguous data.

        Parameters
        ----------
        ext : str
            Extension name suffix, e.g. 'clean', 'nomask', etc.
        im_exts : Iterable(str), None
            Names of image extensions, e.g. 'mask', 'image.pbcor.common', etc.
            If `None` then a default set of extensions is processed. The default
            extensions are "image", "mask", "model", "pb", "psf", "residual".
        """
        if im_exts is None:
            im_exts = (
                    'image', 'image.common', 'image.pbcor.common', 'mask',
                    'model', 'pb', 'psf', 'residual',
            )
        # Validate input.
        im_exts = [im_exts] if isinstance(im_exts, str) else im_exts
        im_exts = [s.lstrip('.') for s in im_exts]
        chunked_imagebases = [c.get_imagebase(ext=ext) for c in self.configs]
        imagebase = self.get_imagebase(ext=ext)
        # Concatenate chunk files for each image extension name.
        for im_ext in im_exts:
            outfile = '{0}.{1}'.format(imagebase, im_ext)
            infiles = [
                    '{0}.{1}'.format(s, im_ext)
                    for s in chunked_imagebases
            ]
            if not os.path.exists(infiles[0]):
                log_post('-- File not found: {0}'.format(infiles[0]), priority='WARN')
                log_post('-- Skipping concatenation for: {0}'.format(im_ext), priority='WARN')
                continue
            log_post(':: Concatenating file "{0}"'.format(outfile))
            if_exists_remove(outfile)
            ia.imageconcat(
                    outfile=outfile,
                    infiles=infiles,
                    reorder=True,
            ).done()
            ia.close()
            ia.done()

    def postprocess(self, ext='clean', use_existing_except=None):
        """
        Parameters
        ----------
        ext : str
        use_existing_except : Iterable(int), None
            If post-processing has already been run, but a few chunks have been
            modified (by manual cleaning for example), then only re-process the
            chunks set by this variable. If left unset (None), then all chunks
            are processed.
        """
        beam = self.get_common_beam(ext=ext)
        assert beam is not None
        imagebase = self.get_imagebase(ext=ext)
        postproc_kwargs = {
                'ext': ext, 'make_fits': False, 'make_hanning': False, 'beam': beam,
        }
        if use_existing_except is None:
            for config in self.configs:
                config.postprocess(**postproc_kwargs)
        else:
            for ix in use_existing_except:
                self.configs[ix].postprocess(**postproc_kwargs)
        self.concat_cubes(ext=ext)
        hanning_smooth_image('{0}.image.common'.format(imagebase))
        export_fits('{0}.image.pbcor.common'.format(imagebase))


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
        do_cont=True, chunked=True):
    """
    Run all pipeline tasks over all setups, spectral windows, and
    uv-weightings.

    Parameters
    ----------
    field : str
    setup : int, None
        Setup number, if `None` run on all setups in sequence.
    weightings : Iterable, default (0.5,)
        List of uv-weightings to use in `tclean`, which may include the string
        "natural" or a number for the briggs robust parameter.
    fullcube : bool, default True
        Image the full spectral window or a small window around the SPWs
        defined rest frequency.
    do_cont : bool, default True
        Whether to image the continuum SPWs or not.
    chunked : bool, default True
        Process frequency intervals of each SPW in serial to reduce memory
        requirements.
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
                if chunked:
                    config.run_pipeline_chunked()
                else:
                    config.run_pipeline()


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
    label : str, None
        SPW label, e.g., "216.113GHz_DCOp". If `None` then apply to all
        Setups and SPWs for target field.
    kind : str
    weighting : str, number
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
# Moment maps
###############################################################################

def make_moments_from_image(imagename, vwin=5, m1_sigma=4, m2_sigma=5,
        overwrite=True):
    """
    Create moment maps for the given image. The cube is masked using the
    associated clean mask.

    Parameters
    ----------
    imagename : str
        CASA image filename ending in '.image'.
    m1_sigma : number
        Multiple of the RMS to threshold the Moment 1 data on.
    m2_sigma : number
        Multiple of the RMS to threshold the Moment 2 data on.
    vwin : number
        Velocity window (half-width) to use for estimating moments over
        relative to the systemic velocity.
        **units**: km/s
    overwrite : bool
    """
    if not os.path.exists(MOMA_DIR):
        os.makedirs(MOMA_DIR)
    # File names
    stem = os.path.splitext(imagename)[0]
    pbname = '{0}.pb'.format(stem)
    maskname = '{0}.mask'.format(stem)
    commonname = '{0}.image.common'.format(stem)
    hannname = '{0}.image.common.hanning'.format(stem)
    basename = os.path.basename(stem)
    outfile = os.path.join(MOMA_DIR, basename)
    # The `ia.moments` task must smooth to a common resolution before
    # computing the moments. To avoid doing this multiple times, process
    # once and store the file.
    if os.path.exists(commonname):
        log_post('-- Using commonbeam image on disk: {0}'.format(commonname))
    else:
        smooth_cube_to_common_beam(imagename)
        copy_pb_mask(commonname, pbname)
    # Calculate Hanning smoothed map for use with higher-order moments.
    if os.path.exists(hannname):
        log_post('-- Using hanning smoothed image on disk: {0}'.format(hannname))
    else:
        hanning_smooth_image(commonname)
    # Create moment maps.
    log_post(':: Calculating moments for: {0}'.format(basename))
    rms = calc_rms_from_image(imagename)
    field = imhead(imagename, mode='get', hdkey='object')
    vsys = ALL_TARGETS[field].vsys
    # Use all emission within the velocity window for m0 and extrema.
    # Moment IDs:
    #    0   moment-0, integrated intensity
    #    8   maximum intensity
    ia.open(commonname)
    imshape = ia.shape()
    region = (
            'box[[0pix,0pix],[{0}pix,{1}pix]], '.format(imshape[0],imshape[1]) +
            'range=[{0:.3f}km/s,{1:.3f}km/s]'.format(vsys-vwin, vsys+vwin)
    )
    ia.moments(
            moments=[0],
            includepix=[0, 1e9],
            region=region,
            outfile=outfile+'.mom0',
            overwrite=True,
    ).done()
    ia.moments(
            moments=[8],
            region=region,
            outfile=outfile+'.max',
            overwrite=True,
    ).done()
    # Use all emission within clean mask for m0 and extrema.
    # Create a T/F mask using an LEL expression on the 1/0 clean mask.
    lel_expr_clean = '"{0}" > 0.5'.format(maskname)
    ia.moments(
            moments=[0],
            mask=lel_expr_clean,
            region=region,
            outfile=outfile+'.mom0_cmask',
            overwrite=True,
    ).done()
    ia.moments(
            moments=[8],
            mask=lel_expr_clean,
            region=region,
            outfile=outfile+'.max_cmask',
            overwrite=True,
    ).done()
    # Apply significance thresholds to m1/m2. Apply hanning smoothing to the
    # spectral axis to help mitigate single-channel noise-spikes.
    # Moment IDs:
    #    1   moment-1, intensity weighted mean velocity
    #    2   moment-2, intensity weighted velocity dispersion
    # Factor of 1.305 reduction in RMS from unsmoothed to Hanning smoothed.
    lel_expr_hann_fmt = '{0} && "{1}" > {{sigma}}*{rms}/1.305'.format(
            lel_expr_clean, hannname, rms=rms,
    )
    ia.moments(
            moments=[1],
            mask=lel_expr_hann_fmt.format(sigma=m1_sigma),
            region=region,
            outfile=outfile+'.mom1',
            overwrite=overwrite,
    ).done()
    ia.moments(
            moments=[2],
            mask=lel_expr_hann_fmt.format(sigma=m2_sigma),
            region=region,
            outfile=outfile+'.mom2',
            overwrite=overwrite,
    ).done()
    ia.done()
    ia.close()
    # Retrieve plane of the primary beam near the systemic velocity
    ia.open(pbname)
    pbcube = ia.getregion(region, dropdeg=True).transpose()
    ia.done()
    ia.close()
    ix_center = pbcube.shape[0] // 2
    pbplane = pbcube[ix_center]
    # Primary beam correct the relevant moment maps
    for ext in ('mom0', 'mom0_cmask', 'max', 'max_cmask'):
        momentname = '{0}.{1}'.format(outfile, ext)
        pbcor_momentname = '{0}.pbcor'.format(momentname)
        impbcor(
                imagename=momentname,
                pbimage=pbplane,
                outfile=pbcor_momentname,
                overwrite=overwrite,
        )
        # Division by NaN appears to create a few "hot pixels" with
        # floating point Inf when dividing by a primary beam plane
        # that has a mask slightly smaller than the mask generated
        # `ia.moments`. Manually replace these Inf values back to NaN.
        ia.open(pbcor_momentname)
        data = ia.getchunk()
        data[np.isinf(data)] = np.nan
        ia.putchunk(data)
        ia.done()
        ia.close()
    # Export the CASA images to FITS files
    export_moment_exts = ('mom0.pbcor', 'mom0_cmask.pbcor', 'max.pbcor',
            'max_cmask.pbcor', 'mom1', 'mom2')
    for ext in export_moment_exts:
        momentname = '{0}.{1}'.format(outfile, ext)
        export_fits(momentname, velocity=True, overwrite=overwrite)


def make_all_moment_maps(field, ext='clean', vwin=5, overwrite=True):
    """
    Generate all moment maps for images with the matching field ID
    name and extension. Moment maps will be written to the directory
    set in ``MOMA_DIR``.

    Parameters
    ----------
    field : str
        Target field ID name
    ext : str
        Image extension name, such as 'clean', 'nomask', etc.
    vwin : number
        Velocity window (half-width) to use for estimating moments over
        relative to the systemic velocity.
        **units**: km/s
    overwrite : bool, default True
        Overwrite moment maps files if they exist.
    """
    image_paths = glob('images/{0}/{0}_*_{1}.image'.format(field, ext))
    image_paths.sort()
    for path in image_paths:
        make_moments_from_image(path, vwin=vwin, overwrite=overwrite)


###############################################################################
# Quality assurance plots
###############################################################################

def savefig(filen, dpi=300, plot_exts=('png', 'pdf'), close=True):
    """
    Save the figure instance to a file.

    Parameters
    ----------
    filen : str
        File name to write to, excluding the extension (e.g., ".pdf").
    dpi : number
    plot_exts : Iterable
        Extensions to create plots for, e.g., "pdf", "svg", "png", "eps", "ps".
    close : bool
        Close figure instance when finished plotting.
    """
    if not os.path.exists(PLOT_DIR):
        os.makedirs(PLOT_DIR)
    fig = plt.gcf()
    for ext in plot_exts:
        # NOTE The GTKAgg backend hits a 16-bit signed integer overflow limit
        # at figure dimensions of approximately 109.23 inches in width or
        # height. This only affects the printing of PNG files, as PDF, PS, and
        # SVG files use different backends that are not limited in this way.
        if ext == 'png' and fig.get_figheight() > 108:
            log_post('-- Skipping PNG, image too large for backend', priority='WARN')
            continue
        plot_filen = os.path.join(PLOT_DIR, '{0}.{1}'.format(filen, ext))
        plt.savefig(plot_filen, dpi=dpi)
    log_post('-- figure saved for: {0}'.format(plot_filen))
    if close:
        plt.close('all')


class CubeSet(object):
    def __init__(self, path, sigma=6):
        """
        Parameters
        ----------
        path : str
            Full image path, including extension.

        sigma : number
            Significance threshold used to select planes/channels to
            retrieve.
        """
        self.path = path
        self.sigma = sigma
        # Setup paths and calculate noise
        stem = os.path.splitext(path)[0]
        self.stem = stem
        self.basename = os.path.basename(stem)
        self.rms = calc_rms_from_image(path)
        # Calculate coordinate offset positions for ticks.
        self.header = imhead(path, mode='list')
        self.pix_scale = abs(np.rad2deg(self.header['cdelt1']) * 60**2)  # arcsec
        # Image size properties
        shape = self.header['shape']
        assert shape[0] == shape[1]  # is square image
        self.shape = shape
        self.pix_width = shape[0]
        self.nchan = shape[3]
        # Identify channels with significant emisssion.
        good_chan = self.get_good_channels()
        self.good_chan = good_chan
        self.ngood = good_chan.size

    @staticmethod
    def get_chunk(imagename):
        ia.open(imagename)
        # Drop (degenerate) stokes axis to reduce dimensions to N=3.
        # Save memory by reducing to 32-bit float. These arrays are only
        # for plotting purposes anyway, so they do not need to preserve
        # full precision.
        chunk = ia.getchunk(dropdeg=True).astype('float32', copy=False)
        ia.done()
        ia.close()
        return chunk.transpose()

    def get_plane_from_image(self, filen, ix):
        """
        Parameters
        ----------
        filen : str
            CASA image path, e.g. "images/CB68/CB68_<...>_clean.mask"
        ix : int
            Channel index number (zero indexed).

        Returns
        -------
        ndarray
            2D image plane of the selected channel.
        """
        try:
            ia.open(filen)
            single_plane = ia.collapse(function='sum', axes=[3], chans=str(ix))
            data = single_plane.getchunk(dropdeg=True)
            single_plane.done()
        except RuntimeError:
            raise
        finally:
            ia.done()
            ia.close()
        # Transpose converts the native Fortran ordering to C ordering.
        return data.astype('float32').transpose()

    def get_image_planes(self, ix):
        stem = self.stem
        im = self.get_plane_from_image(stem+'.image', ix)
        rs = self.get_plane_from_image(stem+'.residual', ix)
        ma = self.get_plane_from_image(stem+'.mask', ix)
        pb = self.get_plane_from_image(stem+'.pb', ix)
        # Mask the image and residual files with NaNs outside of the primary
        # beam limit.
        pbmask = pb < PBLIMIT
        im[pbmask] = np.nan
        rs[pbmask] = np.nan
        return im, rs, ma, pb

    def iter_planes(self):
        for ix in self.good_chan:
            planes = self.get_image_planes(ix)
            yield planes, ix

    def get_good_channels(self):
        stem = self.stem
        threshold = self.sigma * self.rms
        log_post('-- Identifying channels with significant emission')
        imagename = '{0}.image'.format(stem)
        try:
            ia.open(imagename)
            masked_image = ia.imagecalc(
                    pixels='iif("{0}">{1},1,0)'.format(imagename, threshold),
            )
            summed_image = masked_image.collapse(function='sum', axes=[0,1,2])
            data = summed_image.getchunk().squeeze()
            masked_image.done()
            summed_image.done()
        except RuntimeError:
            raise
        finally:
            ia.done()
            ia.close()
        good_chan = np.argwhere(data > 0).flatten()
        return good_chan

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
    Generate Quality Assurance plots by visualizing each channel where
    significant emission occurs in the image cube of interest (i.e.,
    `.image` or `.residual`).

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
        cmap = plt.cm.afmhot
        mask_cont_color = 'cyan'
        vmin, vmax = -3 * cset.rms, 10 * cset.rms
    elif kind == 'residual':
        cmap = plt.cm.RdBu_r
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
            sharey=True, figsize=figsize, dpi=300)
    for ax, (planes, chan_ix) in zip(axes.flat, cset.iter_planes()):
        image, residual, mask, pbeam = planes
        data = image if kind == 'image' else residual
        # Show colorscale of image data
        ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap, origin='lower')
        # Show HPBW contour for primary beam data
        ax.contour(pbeam, levels=[0.5], colors='0.5',
                linestyles='dashed', linewidths=0.4)
        # Show contours for high-SNR emission
        if kind == 'image':
            high_snr_levels = list(cset.rms * np.array([10, 20, 40, 80]))
            ax.contourf(data, levels=high_snr_levels, cmap=plt.cm.rainbow)
        # Show contour for the clean mask
        if mask.any():
            ax.contour(mask, level=[0.5], colors=mask_cont_color,
                    linestyles='solid', linewidths=0.2)
        # show nice channel label in the top-right corner
        text = ax.annotate(str(chan_ix), (0.83, 0.91),
                xycoords='axes fraction', fontsize='xx-small')
        text.set_path_effects([
                path_effects.withStroke(linewidth=2, foreground='white')])
        # show axis ticks as relative offsets in fixed arcsec interval
        ax.set_xticks(tick_pos)
        ax.set_yticks(tick_pos)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        del image, residual, mask, pbeam
    # zoom the window to crop out some of the PB NaNs around the edges
    ax.set_xlim(0.12*cset.pix_width, 0.88*cset.pix_width)
    ax.set_ylim(0.12*cset.pix_width, 0.88*cset.pix_width)
    plt.tight_layout(pad=0.8)
    # Hide unused plots
    for ax in axes.flat[-nempty:]:
        ax.set_visible(False)
    outfilen_ext = '{0}_{1}'.format(outfilen, kind)
    savefig(outfilen_ext)


def make_qa_plots_from_image(path, plot_sigma=None, overwrite=True):
    """
    Create a single set of Quality Assurance plots. Plots will
    be written to the directory specified in ``PLOT_DIR``.

    Parameters
    ----------
    path : str
        Full path to imagename, including the ending ".image".
    plot_sigma : number, None
        Significance threshold to select a channel for plotting. If None, use
        the default set by `CubeSet`.
    overwrite : bool
        Overwrite plot files if they exist.
    """
    # To avoid overwriting files, check if the image PDF file exists
    # and skip plotting if it exists.
    stem = os.path.splitext(path)[0]
    basename = os.path.basename(stem)
    pdf_path = '{0}{1}_qa_plot_image.pdf'.format(PLOT_DIR, basename)
    if not overwrite and os.path.exists(pdf_path):
        log_post('-- File exists, passing: {0}'.format(pdf_path))
        return
    # Read in cube data and make plots
    if plot_sigma is None:
        cset = CubeSet(path)
    else:
        cset = CubeSet(path, sigma=plot_sigma)
    outfilen = '{0}_qa_plot'.format(cset.basename)
    make_qa_plot(cset, kind='image', outfilen=outfilen)
    make_qa_plot(cset, kind='residual', outfilen=outfilen)
    # Delete the CubeSet and trigger a garbage collection to free memory.
    del cset
    gc.collect()


def make_all_qa_plots(field, ext='clean', ignore_chunks=True, overwrite=True):
    """
    Generate all Quality Assurance plots for all image cubes matching
    the given field ID name and extension. Plots will be written to
    the directory set in ``PLOT_DIR``.

    Parameters
    ----------
    field : str
        Target field ID name
    ext : str
        Image extension name, such as 'clean', 'nomask', etc.
    ignore_chunks : bool
        If 'True' ignore chunk image files.
    overwrite : bool, default True
        Overwrite plot files if they exist. Setting to False will avoid the
        potentially large run-time cost of reading cubes into memory to re-make
        existing plots.
    """
    image_paths = glob('{0}{1}/{1}_*_{2}.image'.format(IMAG_DIR, field, ext))
    image_paths.sort()
    # Typical image file path if chunked:
    #   images/CB68/CB68_244.936GHz_CS_chunk3_joint_0.5_clean.image
    pattern = re.compile(r'{0}_.+_.+_chunk[\d]+_'.format(field))
    for path in image_paths:
        if ignore_chunks and pattern.search(path) is not None:
            continue
        make_qa_plots_from_image(path, overwrite=overwrite)
        gc.collect()  # free memory, just in case...


###############################################################################
# User defined functions
###############################################################################

if __name__ == '__main__':
    # NOTE statements placed in this block will be executed on `execfile` in
    #      addition to when called by qsub for batch jobs.
    pass


