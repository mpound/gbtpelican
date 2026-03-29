
import numpy as np
from dysh.fits.gbtfitsload import GBTFITSLoad

sdf=GBTFITSLoad('AGBT25B_386_01/AGBT25B_386_01.raw.vegas')
a=sdf.get_summary()
scans=a[a['OBJECT']=='VANE']['SCAN']
print('VANE scans:',scans)

tsys = np.zeros(16)
for s in scans:
    print(s)
    for b in range(16):
        tsys[b] = sdf.vanecal(scan=s, fdnum=b)
        print(s,b,tsys[b])
    tsys_out = [f"{x:.1f}" for x in tsys]
    print(s, tsys_out)
    print(s, tsys)
