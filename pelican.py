#!/usr/bin/env python
from dysh.fits import GBTFITSLoad
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
        
def pointing_errors():    
    cr2 = si[si['SCAN'].isin([30]) & si['FDNUM'].isin([8]) & si['IFNUM'].isin([0]) & si['PLNUM'].isin([0]) & si['SIG'].isin(['T'])]
    xm=np.nanmean(cr2['CRVAL2'])
    ym=np.nanmean(cr2['CRVAL3'])
    print(f"MEAN RA={xm} DEC={ym}")
    x=(cr2['CRVAL2']-xm).to_numpy()
    y=(cr2['CRVAL3']-ym).to_numpy()
    xyplot(x,y, 'pointing errors scan 30', "arcsec")

def ralongmapview():
    #ralongmap=[31,32,35,36,37,38,39,40,41,51,52,53,54,55,56,57,58]
    cr2 = si[si['PROC'].isin(['RALongMap']) & si['FDNUM'].isin([8]) & si['IFNUM'].isin([0]) & si['PLNUM'].isin([0]) & si['SIG'].isin(['T'])]
    x=(cr2['CRVAL2']-pelicancenter.ra.degree).to_numpy()
    y=(cr2['CRVAL3']-pelicancenter.dec.degree).to_numpy()
    xyplot(x,y, 'RALONGMAP', "arcmin")

def getfs(scan=30):
    sb = sdf.getfs(scan=scan,ifnum=0,plnum=0,fdnum=0)
    return sb

# doesn't work
# do the split with gbtgridder instead
#def split_sb(scanblock):
#    scanblocks = {} 
#    scanblocks["HCO+"] = scanblock.copy()
#    scanblocks["HCN"] = scanblock.copy()
##    for k in scanblocks: 
#        #rf= restfreq[k]    
#        _cc = scanblocks[k][0]._calibrated.copy()
#        scanblocks[k][0]._calibrated = _cc[:,channels[k]]
#    return scanblocks

def specplot(scan=30):
    sb = getfs(scan)
    spec = sb.timeaverage()
    p = spec.plot()
    p.show_catalog_lines(cat='gbtlines',chemical_name="Formylium")
    p.show_catalog_lines(cat='gbtlines',chemical_name="Hydrogen Cyanide")

def load(file,loadAll = True):
    _sdf = GBTFITSLoad(file)
    if loadAll:
        print(f"loading all of {file}")
        _sdf.load_all()
    return _sdf

#def createScanBlock(line:str, sdf:GBTFITSLoad, scans:list, channel:list ) -> ScanBlock:
#    sb = sdf.getfs(scan=scans,ifnum=0,plnum=0,fdnum=0,channel=channel)
#    for i in range(0,15):
#        sb[k].append(sdf.getfs(scan=scans,ifnum=0,plnum=0,fdnum=i,channel=channels[k]))
class BaselineAvg(NamedTuple):
    timeaverage:list = []
    model: Model = None
    line: str = None
    
def baseline(scanblock:list,subtract:bool, restvalue:Quantity, line: str, plot=False):#  -> BaselineAvg:
    # scanblock i contains all integrations for all scans of feed i
    # We want to compute the average spectrum for all integrations so we can find the good baseline.
    
    tavg=[]

    for i in range(len(scanblock)):
        print(f"Doing baseline {i} of len(scanblock)")
        tavg.append(scanblock[i].timeaverage(use_wcs=False))
        print(f"fdnum={i}, exp={tavg[i].meta["EXPOSURE"]}")
        tavg[i].rest_value=restvalue
        tavg[i].baseline(degree=2,exclude=[60,100],remove=True)
        #lsrk=tavg[i].with_frame("LSRK")
        # bad channel at 30.25 km/s
        #lsrk.baseline(degree=2,include=[(-20*kms,-6*kms),(3*kms,20*kms),(30*kms,31*kms)],remove=True)
        if plot:
            #lsrk.plot(ymin=-1,ymax=1,yaxis_unit="K")
            tavg[i].with_frame('LSRK').plot(yaxis_unit="mK",xaxis_unit='km/s')
        if subtract:
            scanblock[i].undo_baseline()
            scanblock[i].subtract_baseline(tavg[i].baseline_model,tol=1E5)
        #return BaselineAvg(timeaverage=tavg,  model=lsrk.baseline_model, line=line)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    parser.add_argument("--feeds",type=int,nargs='+', action='store',help="list feeds to reduce, e.g. --feeds 0 1 2 3 ", required=True)
    parser.add_argument("--line",type=str, action='store',help="which line to reduce. One of HCO+, HCN (case insensisitive) ", required=True)
    parser.add_argument("--file", type=str, default=None, action='store', help='output file for final scanblock')
    parser.add_argument('--baseline',default=False,action='store_true',help='remove baseline')
    parser.add_argument('--plot',default=False,action='store_true',help='plot removed baseline')
    parser.add_argument('--verbose', '-v', default=False, action='store_true', help='verbose mode')
    parser.add_argument('--write', '-w', default=True, action='store_false', help='write to pelican_line_feed.sdfits')
    args = parser.parse_args()
    
    
    lines = ["HCO+","HCN"]
    args.line = args.line.upper()
    if args.line not in lines:
        raise ValueError(f"line must be one of {lines}")
    channels = {}
    startchan = {}
    final_sb={}
    restfreq = {}
    restfreq["HCO+"] = 89.188518*u.GHz
    restfreq["HCN"] = 88.6318473*u.GHz

    # use a total width of 0.01 GHz centered on rest frequencies
    startchan["HCO+"] = 3800
    startchan["HCN"] = 12311
    nchan = 151
    final_sb["HCO+"] = ScanBlock()
    final_sb["HCN"] = ScanBlock()
    kms=u.km/u.s
    pelicancenter=SkyCoord('20:51:08.07 +44:26:35.3',frame='fk5',unit=(u.hr,u.degree))
    
    for k in restfreq:
        channels[k] = [startchan[k],startchan[k]+nchan]

    fnm = "/bigdisk/data/gbt/AGBT25B_386_01/AGBT25B_386_01.raw.vegas/"
    print(f"Loading {fnm}")
    sdf = load(fnm,loadAll=False)

    si = sdf._sdf[0]._index
    scans=[31,32,35,36,37,38,39,40,41,51,52,53,54,55,56,57,58]
    #scans=[30]
    maxfeed=15
    feeds = np.arange(0,maxfeed)
    #@TODO change sb to ScanBlock and use sb.extend to concat the feeds

    doplot=args.plot
    dobase=args.baseline
    lsrk = {"HCO+": [], "HCN":[]}
    print(f"Doing feed={args.feeds}")
    feedcount = 0
    for k in [args.line]:
        for i in args.feeds:
            print(f"Doing {k} feed {i} {scans=} channels={channels[k]}")
            final_sb[k].append(sdf.getfs(scan=scans,ifnum=0,plnum=0,fdnum=i,channel=channels[k]))
        if dobase:
            baseline(final_sb[k],subtract=True,restvalue=restfreq[k],line=k,plot=doplot)
        for i in args.feeds:          
            _x = final_sb[k][feedcount].timeaverage()
            _x.rest_value=restfreq[k]
            lsrk[k].append(_x.with_frame('LSRK'))
            feedcount= feedcount+1
        if args.write:
            count = 0
            for i in args.feeds:
                final_sb[k][count].write(f"pelican_{args.line}_{i}.sdfits",flags=True)
                count = count+ 1
        
    sball = {"HCN": ScanBlock(), "HCO+": ScanBlock()}
    if len(args.feeds) > 1:
        for i in range(len(final_sb[args.line])):
            sball[args.line].extend(final_sb[args.line][i])
        final_avg = sball[args.line].timeaverage()
    if args.plot:
        sball[args.line].timeaverage().with_frame('LSRK').plot(xaxis_unit='km/s')
    if args.write:
        sball[args.line].write('pelican_{args.line}_all.sdfits',flags=True)
        
    #final_sb.plot(vmin=-0.005,vmax=0.005)
        #tavg[i].with_frame("LSRK").plot(xaxis_unit="km/s",xmin=-15,xmax=15)
        #tavg[i].rest_value=restfreq["HCN"]
        #tavg[i].with_frame("LSRK").plot(xaxis_unit="km/s",xmin=-15,xmax=15)
