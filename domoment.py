#!/usr/bin/env python
"""Command-line wrapper around maskmoment.maskmoment (which has no CLI of its own)."""
import argparse
import os
import sys

# maskmoment is an add-on git clone, not installed; make it importable.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "maskmoment"))
from maskmoment import maskmoment  # noqa: E402

parser = argparse.ArgumentParser(description="Make masked moment maps of a cube with maskmoment.")
parser.add_argument("cube", help="FITS image cube (K or Jy/beam)")
parser.add_argument("--gain_fits", default=None, help="gain cube")
parser.add_argument("--rms_fits", default=None, help="noise cube")
parser.add_argument("--rms", type=float, default=None, help="global rms noise, in cube units")
parser.add_argument("--outdir", default="", help="output directory. Default: where the cube is")
parser.add_argument("--outname", default=None, help="output basename. Default: root of cube name")
parser.add_argument("--snr_hi", type=float, default=4, help="high significance threshold (default 4)")
parser.add_argument("--snr_lo", type=float, default=2, help="low significance threshold (default 2)")
parser.add_argument("--snr_hi_minch", type=int, default=1)
parser.add_argument("--snr_lo_minch", type=int, default=1)
parser.add_argument("--min_tot_ch", type=int, default=2)
parser.add_argument("--minbeam", type=float, default=1)
parser.add_argument("--nguard", type=int, nargs=2, default=[0, 0], help="sky pixels, channels")
parser.add_argument("--edgech", type=int, default=5, help="edge channels used for rms")
parser.add_argument("--fwhm", type=float, default=None, help="spatial smoothing FWHM (arcsec)")
parser.add_argument("--vsm", type=float, default=None, help="spectral smoothing width (channels)")
parser.add_argument("--vsm_type", default="gauss", choices=["boxcar", "gauss", "gaussfinal"])
parser.add_argument("--vmin", type=float, default=None, help="min channel for moments")
parser.add_argument("--vmax", type=float, default=None, help="max channel for moments")
parser.add_argument("--output_peak", action="store_true", help="also output peak, vpeak, effective width")
args = parser.parse_args()

kwargs = vars(args)
cube = kwargs.pop("cube")
maskmoment(cube, **kwargs)
