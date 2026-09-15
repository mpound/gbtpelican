#!/usr/bin/env python
from dysh.fits import GBTFITSLoad
#from dysh.line import SpectralLineSearch
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

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
    
restfreq = {}
restfreq["HCO+"] = 89.188518*u.GHz
restfreq["HCN"] = 88.6318473*u.GHz
channels = {}
startchan = {}
nchan = 151
# use a total width of 0.01 GHz centered on rest frequencies
startchan["HCO+"] = 3917
startchan["HCN"] = 12361
for k in restfreq:
    channels[k] = [startchan[k],startchan[k]+nchan]

fnm = "/bigdisk/data/gbt/AGBT25B_386_01/AGBT25B_386_01.raw.vegas/"
sdf = load(fnm,loadAll=False)
pelicancenter=SkyCoord('20:51:08.07 +44:26:35.3',frame='fk5',unit=(u.hr,u.degree))
si = sdf._sdf[0]._index
scans=[31,32,35,36,37,38,39,40,41,51,52,53,54,55,56,57,58]
sb = []
tavg=[]
do="HCO+"
kms=u.km/u.s
for i in range(0,15):
    sb.append(sdf.getfs(scan=scans,ifnum=0,plnum=0,fdnum=i))
for i in range(0,15):
    #tavg.append(sb[i].timeaverage(use_wcs=False))
    print(f"fdnum={i}, exp={tavg[i].meta["EXPOSURE"]}")
    tavg[i].rest_value=restfreq[do]
    lsrk=tavg[i].with_frame("LSRK")
    lsrk.baseline(degree=2,include=[(-20*kms,-6*kms),(3*kms,20*kms)],remove=True)
    #lsrk.plot(ymin=0.002,ymax=0.003,yaxis_unit="K")
    sb[i].subtract_baseline(lsrk.baseline_model,tol=1E5)
    
sb[6].plot(vmin=-0.005,vmax=0.005)
    #tavg[i].with_frame("LSRK").plot(xaxis_unit="km/s",xmin=-15,xmax=15)
    #tavg[i].rest_value=restfreq["HCN"]
    #tavg[i].with_frame("LSRK").plot(xaxis_unit="km/s",xmin=-15,xmax=15)
