#!/usr/bin/env python
from dysh.fits import GBTFITSLoad, GBTOffline
from dysh.spectra.scan import ScanBlock
#from dysh.line import SpectralLineSearch
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import Model
from typing import NamedTuple
from astropy.units.quantity import Quantity
import argparse
import sys
#import os
from astropy.table import Table

def xyplot(xdeg:np.ndarray,ydeg:np.ndarray,title:str,unit:str|u.Unit):
    x=(xdeg*u.degree).to(unit)
    y=(ydeg*u.degree).to(unit)
    plt.title(title)
    plt.xlabel(fr'$\Delta\alpha$ ({unit})')
    plt.ylabel(fr'$\Delta\delta$ ({unit})')
    plt.gca().invert_xaxis()
    plt.plot(x,y)
    plt.show()

def plotbeams(insdf,fdnum=None,unit:str|u.Unit="arcmin"):
    df=insdf[["SCAN","FDNUM","CRVAL2","CRVAL3","PLNUM","IFNUM","PROC"]]
    df = df[df["PROC"].isin(["RALongMap"])]
    if fdnum is not None: 
        df = df[df["FDNUM"].isin(fdnum)]
        title = f"Beam(s) {fdnum}"
    else:
        title = "Beams"
    x=(df['CRVAL2']-pelicancenter.ra.degree).to_numpy()*u.degree
    y=(df['CRVAL3']-pelicancenter.dec.degree).to_numpy()*u.degree
    x=x.to(unit)
    y=y.to(unit)

    plt.xlabel(fr'$\Delta\alpha$ ({unit})')
    plt.ylabel(fr'$\Delta\delta$ ({unit})')
    plt.scatter(x,y,marker='o',s=1)
    plt.title(title)
    plt.show()
        
def pointing_errors(si):    
    cr2 = si[si['SCAN'].isin([30]) & si['FDNUM'].isin([8]) & si['IFNUM'].isin([0]) & si['PLNUM'].isin([0]) & si['SIG'].isin(['T'])]
    xm=np.nanmean(cr2['CRVAL2'])
    ym=np.nanmean(cr2['CRVAL3'])
    print(f"MEAN RA={xm} DEC={ym}")
    x=(cr2['CRVAL2']-xm).to_numpy()
    y=(cr2['CRVAL3']-ym).to_numpy()
    xyplot(x,y, 'pointing errors scan 30', "arcsec")

def ralongmapview(si):
    #ralongmap=[31,32,35,36,37,38,39,40,41,51,52,53,54,55,56,57,58]
    cr2 = si[si['PROC'].isin(['RALongMap']) & si['FDNUM'].isin([8]) & si['IFNUM'].isin([0]) & si['PLNUM'].isin([0]) & si['SIG'].isin(['T'])]
    x=(cr2['CRVAL2']-pelicancenter.ra.degree).to_numpy()
    y=(cr2['CRVAL3']-pelicancenter.dec.degree).to_numpy()
    xyplot(x,y, 'RALONGMAP', "arcmin")

def getfs(sdf, scan):
    sb = sdf.getfs(scan=scan,ifnum=0,plnum=0,fdnum=0)
    return sb

def specplot(scan=30):
    sb = getfs(scan)
    spec = sb.timeaverage()
    p = spec.plot()
    p.show_catalog_lines(cat='gbtlines',chemical_name="Formylium")
    p.show_catalog_lines(cat='gbtlines',chemical_name="Hydrogen Cyanide")
def loadsession(file,loadAll = False):
    _sdf = GBTOffline(file)
    print(f"{file} loaded")
    if loadAll:
        print(f"loading all of {file}")
        _sdf.load_all()
    return _sdf

def load(file,loadAll = True):
    _sdf = GBTFITSLoad(file)
    if loadAll:
        print(f"loading all of {file}")
        _sdf.load_all()
    return _sdf

class BaselineAvg(NamedTuple):
    timeaverage:list = []
    model: Model = None
    line: str = None
    
def baseline_one(sb, degree:int, subtract:bool, restvalue:Quantity, exclude=[[0,5],[60,100]]):
    """Fit a baseline to the time average of one feed's ScanBlock, optionally
    subtracting it from the ScanBlock. Does no plotting.
    Returns the baseline-fitted time average."""
    tavg = sb.timeaverage(use_wcs=False)
    tavg.rest_value=restvalue
    tavg.baseline(degree=degree,exclude=exclude,remove=True)
    if subtract:
        sb.undo_baseline()
        sb.subtract_baseline(tavg.baseline_model,tol=1E5)
    return tavg

def baseline(scanblock:list,degree:int, subtract:bool, restvalue:Quantity, exclude=[[0,5],[60,100]], plot=False):
    # scanblock i contains all integrations for all scans of feed i
    # We want to compute the average spectrum for all integrations so we can find the good baseline.
    tavg=[]
    for i in range(len(scanblock)):
        print(f"Doing baseline {i} of {len(scanblock)}")
        tavg.append(baseline_one(scanblock[i],degree,subtract,restvalue,exclude))
        print(f'fdnum={i}, exp={tavg[i].meta["EXPOSURE"]}')
        if plot:
            tavg[i].with_frame('LSRK').plot(ymin=-2,ymax=2,yaxis_unit="mK",xaxis_unit='km/s')
    return tavg

def reduce_feed(sdf, fdnum:int, scans, units:str, zenith_opacity, channels, baseorder:int, restvalue:Quantity):
    """Calibrate (and optionally baseline-subtract) one feed. Returns (ScanBlock, tavg or None)."""
    print(f"Doing feed {fdnum} scans={scans} channels={channels}")
    sb = sdf.getfs(scan=scans,ifnum=0,plnum=0,fdnum=fdnum, units=units, zenith_opacity=zenith_opacity,
                   channel=channels,shift_method='interpolate')
    tavg = None
    if baseorder >= 0:
        tavg = baseline_one(sb, baseorder, True, restvalue)
        print(f'fdnum={fdnum}, exp={tavg.meta["EXPOSURE"]}')
    return sb, tavg

class LineParams:
    """Parameters for each Line map/reduction"""
    def __init__(self, line:str, restfreq:Quantity, sdf:GBTFITSLoad, startchan:int = 1, nchan: int = 16384,scans=[]):
        self.line = line.upper()
        self.restfreq = restfreq
        self.startchan = startchan
        self.nchan = nchan
        self.channels = (startchan,startchan+nchan)
        self.final_sb = ScanBlock()
        self.scans = scans     
    @property
    def endchan(self):
        """end channel"""
        return self.channels[-1]

    def concat(self):
        """Concatenate all ScanBlocks"""
        s = ScanBlock()
        for f in self.final_sb:
            s.extend(f)
        return s

    def write(self, file, **kwargs ):
        self.concat().write(file, **kwargs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    #parser.add_argument('--baseline',default=False,action='store_true',help='remove baseline')
    # @TODO allow negatives to be "all but -M -N feeds"
    parser.add_argument("--feeds",type=int,nargs='+', action='store',help="list feeds to reduce, e.g. --feeds 0 1 2 3 ", required=True)  
    parser.add_argument("--input", type=str, default=None, action='store', help='tabular input file with Project, Line, RestFreq, scans etc.')
    parser.add_argument('--flux',  action="store_true", help='Use Flux(Jy) instead of Ta(K)')
    parser.add_argument("--line",type=str, action='store',help="which line to reduce. One of HCO+, HCN (case insensisitive) ", required=True)
    parser.add_argument("--nchan",type=int, default=151,help="how many channels to keep")
    parser.add_argument('--plot',default=False,action='store_true',help='plot spectra')
    parser.add_argument('--verbose', '-v', default=False, action='store_true', help='verbose mode')
    parser.add_argument('--session', default="all", type=str, action='store', help='which session to reduce. Default: all sessions input file')
    parser.add_argument('--split', default=False, action='store_true', help='If true, split sdfits output files into separate feeds')
    parser.add_argument('--output', '-o',action='store', help='File to write to. Default: pelican_{line}_{session}.sdfits')
    args = parser.parse_args()
    

    # Constant parameters
    kms=u.km/u.s
    pelicancenter=SkyCoord('20:51:08.07 +44:26:35.3',frame='fk5',unit=(u.hr,u.degree))
    
    table = Table.read(args.input, format='ascii', 
                       converters={'Session':str, 'RestFreq(GHz)':np.double})
    lines=set(table['Line'])
    if args.line.upper() not in lines:
        raise ValueError(f"line must be one of {lines}")

    for i in range(len(table)):
        line = table['Line'][i]
        # only do the requested line
        if(line.upper() != args.line.upper()):
            continue
        session = table['Session'][i]
        if (args.session != "all") and args.session != session:
            if args.verbose:
                print(f"Skipping session {session}")
            continue
        filesession = f"{table['Project'][i]}_{session}"
        sdf = loadsession(filesession)
        startchan = table['startchan'][i]
        zenith_opacity = table['tau'][i]
        nchan = table['nchan'][i]
        s = table['scans'][i]
        scans = [int(x) for x in s.split(',')]
        baseorder = table['baseline'][i]
        rf = table['RestFreq(GHz)'][i]*u.GHz
        lp = LineParams(line,rf,sdf, startchan, nchan, scans=scans)
        doplot=args.plot
        dobase  = (baseorder >=0 )
        lsrk = []
        print(f"Doing feed={args.feeds}")
        feedcount = 0
        if args.flux:
            units="flux"
        else:
            units="ta"
        feedargs = (scans, units, zenith_opacity, lp.channels, baseorder, lp.restfreq)
        results = [reduce_feed(sdf, fd, *feedargs) for fd in args.feeds]
        for _sb, _ in results:
            lp.final_sb.append(_sb)
        if doplot:
            for _, tavg in results:
                if tavg is not None:
                    tavg.with_frame('LSRK').plot(ymin=-2,ymax=2,yaxis_unit="mK",xaxis_unit='km/s')
        if doplot:
            scanplots = [x.plot(vmin=-0.005,vmax=0.005) for x in lp.final_sb]
        if False:
            for fdnum in args.feeds:          
                _x = lp.final_sb[feedcount].timeaverage()
                _x.rest_value=lp.restfreq
                lsrk.append(_x.with_frame('LSRK'))
                feedcount= feedcount+1
        if args.output is None:
            lp.write(f"pelican_{args.line}_Session{table['Session'][i]}.sdfits", flags=False, overwrite=True)
        else:
            lp.write(args.output, flags=False, overwrite=True)

# ./gbtgridder/src/gbtgridder --scans 60,61,62,63,64,65,66,77,78 --restfreq 97980.96 --mapcenter 312.783625 44.44313889 --output pelican_CS_05 pelican_CS_Session05.sdfits
