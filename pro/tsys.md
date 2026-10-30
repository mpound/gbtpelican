# Goal
Identify poorly performing feeds/beams in data by examining their system temperatures compared to each other.

## Context

The Argus Receiver on the GBT has 16 beams. These are referred
to as feed number 0 through 15 or beam 1 through 16.   The files
dobeams_AGBT25B_386\*.log contain measurements of system temperatures
for sessions numbered 01,03,04,05,07. There are several measurements over the course of a night indicated by Scan number. In many cases 1 or more beams has
a high system temperature compared to others measured at the same time
indicating it is performing poorly.  A typical good system temperature is less than about 300 K.

Note there may be cases where all beams appeared bad, with extremely
high or extremely negative values, or values close to zero. These were
bad measurements and should be ignored.

## Desired output
A table indicating for each sessions what beams appeared underperforming and for what fraction of the scans in that session it was underperforming. Include also columns with feed number and average Tsys for poorly performing scans.

## Desired implementation

1. Modify ../tsys.py to parse the output log file if args.beams is True.
2. Add an option to tsys.py to parse existing log files.

## Implemented options (tsys.py)

| Option | Default | Meaning |
|---|---|---|
| `--beams` / `--no-beams` | on | Run `vanecal` in GBTIDL for per-beam Tsys (writes `pro/dobeams_<project>.log`), then print the bad-beam table. |
| `--parse-only` | off | Skip all GBTIDL work; parse existing `pro/dobeams_<project>.log` files and print the table. Works off-site. |
| `--bad-factor K` | 1.3 | Flag a beam as underperforming if its Tsys > K * the median of all 16 beams in that scan. Beams with Tsys < `--tsys-min` are also flagged. |
| `--tsys-min` | 50 K | Reject a whole scan if its median beam Tsys is below this. |
| `--tsys-max` | 1000 K | Reject a whole scan if its median beam Tsys is above this. Scans with any non-finite value are also rejected. |

Example, using existing logs:

    ./tsys.py --parse-only AGBT25B_386_01 AGBT25B_386_03 AGBT25B_386_04 AGBT25B_386_05 AGBT25B_386_07

Incomplete scans (fewer than 16 beams in the log) are ignored.

## Output table

One block per session, one row per underperforming feed, sorted by fraction of scans affected:

| Column | Meaning |
|---|---|
| Session | Project/session name. |
| Feed, Beam | Feed number (0-15) and beam number (1-16). |
| Nbad/Nscan | Scans in which the feed was flagged / scans used (after rejection). |
| Frac | Nbad/Nscan. |
| med(Tsys)_bad | Median Tsys of the feed over its flagged scans. |
| med(Tsys)_good | Median, over the same flagged scans, of the per-scan median of the unflagged beams. A baseline for comparison with med(Tsys)_bad. |

The line above each block gives the number of scans used and rejected.


As Marc if you have questions. 
