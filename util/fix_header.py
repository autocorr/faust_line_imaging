#!/usr/bin/env python3
"""
Fix Header
==========
Fix FITS header keywords for ALMA archive delivery. Sets ``DATAMIN``,
``DATAMAX``, ``INSTRUME``, and ``OBJECT`` in the primary HDU header,
updating the file in-place.

To update a set of files call this script from the command line or,
alternatively, import this module and loop through the files. For
example, to the fix header files for L1527 image products, call the
following from the command line:

```
$ # For a single file:
$ python <path>/<to>/util/fix_header.py L1527 L1527_262.004GHz_CCH_joint_0.5_clean.max.pbcor.fits
Processing:L1527_262.004GHz_CCH_joint_0.5_clean.max.pbcor.fits
$ # For multiple files within a directory:
$ python <path>/<to>/util/fix_header.py L1527 L1527*.fits
```

Note that this requires ``astropy`` and ``numpy`` to be installed
into the current Python environment.
"""

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def fix_header(fits_path, source_id):
    """
    Fix FITS header keywords for ALMA archive delivery.

    Opens the FITS file in update mode and sets ``DATAMIN``, ``DATAMAX``,
    ``INSTRUME``, and ``OBJECT`` in the primary HDU header.

    Parameters
    ----------
    fits_path : pathlib.Path or str
        Path to the FITS file to update in-place.
    source_id : str
        Source name written to the OBJECT keyword (e.g. 'L1527').
    """
    fits_path = Path(fits_path)
    with fits.open(fits_path, mode='update') as hdul:
        header = hdul[0].header
        data = hdul[0].data
        header['DATAMIN'] = (float(np.nanmin(data)), '')
        header['DATAMAX'] = (float(np.nanmax(data)), '')
        header['INSTRUME'] = ('ALMA', '')
        header['OBJECT'] = (source_id, '')
        hdul.flush()


def main():
    parser = argparse.ArgumentParser(
        description='Fix FITS headers for ALMA archive delivery.')
    parser.add_argument('source_id', help='Source name (e.g. L1527)')
    parser.add_argument('fits_files', nargs='+', type=Path,
                        help='One or more FITS files to update in-place')
    args = parser.parse_args()
    for p in args.fits_files:
        print(f'Processing: {p.name}')
        fix_header(p, args.source_id)


if __name__ == '__main__':
    main()
