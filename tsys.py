#! /usr/bin/env python
#
#      compute <Tsys> for all VANES in a GBT session

"""
   Example Usage:
       ./tsys.py GBT25B_386_01


   Produces pro/AGBT21B_024_06.tsys which looks something like this:

#   SCAN       Tsys_mean       Tsys_rms        Tsys_min        Tsys_max
      17       181.86916       14.323252       158.81550       216.11503
      54       178.42336       13.545066       155.68764       210.57443
      61       180.99135       14.369497       156.98510       214.64585
      98       183.66805       14.172011       159.27379       217.18325
     105       165.61667       16.774871       146.14036       216.92920
     142       169.82475       17.154855       149.82684       222.31952

    Timing:  AGBT21B_024_02 took 38 mins on fourier  (29 VANE's)

    Multiple sessions can be given in the command line arguments

    NOTE: this command can only run at GBT
"""

import os
import sys
import math
import statistics
import argparse

parser = argparse.ArgumentParser(description="Compute Tsys for all VANEs in GBT sessions")
parser.add_argument("projects", metavar="project", nargs="+",
                    help="GBT session(s), e.g. AGBT25B_386_01")
parser.add_argument("--tsys", action=argparse.BooleanOptionalAction, default=True,
                    help="run vanecal2 and write the .tsys file (default: %(default)s)")
parser.add_argument("--beams", action=argparse.BooleanOptionalAction, default=True,
                    help="run vanecal for per-beam Tsys (default: %(default)s)")
parser.add_argument("--vc3", action=argparse.BooleanOptionalAction, default=False,
                    help="run vanecal3 (default: %(default)s)")
parser.add_argument("--freeze", action=argparse.BooleanOptionalAction, default=False,
                    help="issue 'freeze' in the GBTIDL scripts to stop plotter from updating (default: %(default)s)")
parser.add_argument("--parse-only", action="store_true",
                    help="skip GBTIDL; only parse existing pro/dobeams_<project>.log files and print the bad-beam table")
parser.add_argument("--bad-factor", type=float, default=1.3,
                    help="flag a beam if Tsys > factor * scan median (default: %(default)s)")
parser.add_argument("--tsys-min", type=float, default=50.0,
                    help="reject a scan if its median beam Tsys is below this, in K (default: %(default)s)")
parser.add_argument("--tsys-max", type=float, default=1000.0,
                    help="reject a scan if its median beam Tsys is above this, in K (default: %(default)s)")
args = parser.parse_args()

dotsys = args.tsys
dobeams = args.beams
dovc3 = args.vc3
dofreeze = args.freeze


def parse_beams(logfile):
    """Parse a dobeams log into {scan: [(feed, Tsys), ...]}, keeping only complete 16-beam scans."""
    scans = {}
    scan = None
    with open(logfile) as fp:
        for line in fp:
            words = line.split()
            if words[:1] == ['Scan:']:
                scan = int(words[1])
                scans[scan] = []
            elif scan is not None and line.startswith('fdnum, beam'):
                scans[scan].append((int(words[3]), float(words[5])))
    return {s: b for s, b in scans.items() if len(b) == 16}


def beam_table(projects, factor, tmin, tmax):
    """Print, per session, the beams exceeding factor * scan median and how often."""
    print("# Session           Feed Beam  Nbad/Nscan  Frac  med(Tsys)_bad[K]")
    for project in projects:
        logfile = "dobeams_%s.log" % project
        if not os.path.exists(logfile):
            print("# %s: %s not found" % (project, logfile))
            continue
        scans = parse_beams(logfile)
        good = {}
        for scan, beams in scans.items():
            vals = sorted(t for f, t in beams)
            med = (vals[7] + vals[8]) / 2
            if all(math.isfinite(t) for t in vals) and tmin <= med <= tmax:
                good[scan] = (beams, med)
        print("# %s: %d scans used, %d rejected" % (project, len(good), len(scans) - len(good)))
        bad = {}
        for scan, (beams, med) in good.items():
            for feed, t in beams:
                if t > factor * med or t < tmin:
                    bad.setdefault(feed, []).append(t)
        rows = sorted(bad.items(), key=lambda kv: -len(kv[1]))
        for feed, ts in rows:
            print("  %-18s %4d %4d  %4d/%-5d %5.2f  %8.1f" %
                  (project, feed, feed + 1, len(ts), len(good), len(ts) / len(good), statistics.median(ts)))


#                  our GBTIDL procedures live in pro/
print("### Working in pro")
os.chdir('pro')

if args.parse_only:
    beam_table(args.projects, args.bad_factor, args.tsys_min, args.tsys_max)
    sys.exit(0)

for project in args.projects:
    #  @todo rawdata/project must exist, as a check
    print("### Working on %s" % project)


    print("### writing pro/do1.pro to set the summary")
    fp = open('do1.pro','w')
    fp.write('offline,"%s"\n' % project)
    fp.write('summary,"%s.summary"\n' % project)
    fp.write('exit\n')
    fp.close()


    print("### gbtidl do1.pro > do1.log")
    os.system('gbtidl do1.pro > do1.log')


    print("### parsing pro/%s.summary to find the VANE scans" % project)
    vanes = []
    fp = open("%s.summary" % project)
    lines = fp.readlines()
    for line in lines:
        words = line.strip().split()
        if len(words) > 2:
            if words[1] == 'VANE':
                vanes.append(words[0])
    fp.close()
    nvanes = len(vanes)

    print("### writing pro/do2.pro")
    fp = open('do2.pro','w')
    fp.write('offline,"%s"\n' % project)
    if dofreeze:
        fp.write('freeze\n')
    for v in vanes:
        fp.write('vanecal2,%s\n' % v)
    fp.write('exit\n')
    fp.close()


    print("### gbtidl do2.pro > do2.log")
    os.system('gbtidl do2.pro > do2.log')


    if dotsys:
        print("### parsing pro/do2.log to find the tsys lines for the %d VANES" % nvanes)
        fp = open("do2.log")
        lines = fp.readlines()
        grab = False
        tsys = []
        tau0 = ''
        for line in lines:
            words = line.strip().split()
            if len(words) > 2:
                if words[0][:8] == '(zenith)':
                    tau0 = words[3]
                if grab:
                    tsys.append(line.strip() + ' ' + tau0 + '\n')
                    grab = False
                if words[0][:4] == 'Tatm':
                    grab = True
        fp.close()


        print("### writing pro/%s.tsys with %d vanes" % (project,nvanes))
        fp = open("%s.tsys" % project, "w")
        fp.write("#  %s\n" % project)
        #fp.write("#   SCAN       Tsys_mean       Tsys_rms        Tsys_min        Tsys_max\n")
        #          17  307.74   23.20  270.93  357.18  437.21  372.15  257.60  
        fp.write("# SCAN Tsys_mean rms min     max     Tcal    Twarm   Tatm  tau0\n")
        for t in tsys:
            fp.write(t)
        fp.close()

    if dovc3:
        print("### writing pro/do3.pro")
        fp = open('do3.pro','w')
        fp.write('offline,"%s"\n' % project)
        if dofreeze:
            fp.write('freeze\n')
        for v in vanes:
            fp.write('vanecal3,%s\n' % v)
        fp.write('exit\n')
        fp.close()

        outcmd = f"gbtidl do3.pro> do3_{project}.log"
        print(f"### {outcmd}")
        os.system(outcmd)

    if dobeams:
        print("### writing pro/dobeams.pro")
        fp = open('dobeams.pro','w')
        fp.write('offline,"%s"\n' % project)
        if dofreeze:
            fp.write('freeze\n')
        for v in vanes:
            fp.write('vanecal,%s\n' % v)
        fp.write('exit\n')
        fp.close()

        outcmd = f"gbtidl dobeams.pro> dobeams_{project}.log"
        print(f"### {outcmd}")
        os.system(outcmd)

if dobeams:
    beam_table(args.projects, args.bad_factor, args.tsys_min, args.tsys_max)

print("All done!")


"""

Connecting to file: /home/sdfits/AGBT21B_024_02/AGBT21B_024_02.raw.vegas
Scan:    12 (IF:0 FD:9 PL:0)    Tsys:   1.00
(zenith) Opacity(113.52528,59527.096) = 0.3426
AtmTsys(113.52528,59527.096) = 80.4928
Tatm(113.52528,59527.096) = 269.2360
      12       188.04533       15.111458       161.62504       225.54843
Scan:    18 (IF:0 FD:9 PL:0)    Tsys:   1.00
(zenith) Opacity(113.52528,59527.100) = 0.3427
AtmTsys(113.52528,59527.100) = 80.3749
Tatm(113.52528,59527.100) = 269.1996
      18       190.07192       16.442001       161.41116       229.50203
Scan:    55 (IF:0 FD:9 PL:0)    Tsys:   1.00
(zenith) Opacity(113.52527,59527.119) = 0.3430
AtmTsys(113.52527,59527.119) = 80.2163
Tatm(113.52527,59527.119) = 269.0268
      55       191.17041       16.998825       162.73008       233.13705
"""
