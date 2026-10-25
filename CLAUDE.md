# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Scripts for reducing GBT project AGBT25B_386 (Pelican Nebula, ARGUS 16-beam receiver, HCO+ and HCN zero-spacing maps to complement CARMA/BIMA). Not a package: a handful of standalone scripts, a Makefile for GBO/session housekeeping, and IDL (GBTIDL) procedures. No test suite or linter. `test_pelican.ipynb` is an exploratory notebook.

## Main script: `pelican.py`

Reduces one spectral line for selected feeds using `dysh` (`GBTFITSLoad`, `ScanBlock`), writing an SDFITS file.

```
./pelican.py --feeds 0 1 2 --line HCO+ --input input.tab [-o out.sdfits] [--flux] [--plot] [--split]
```

- Per-session reduction parameters live in a whitespace-delimited table (`input.tab`, `#` comments): Project, Date, Session, Line, RestFreq(GHz), startchan, nchan, baseline order, tau (zenith opacity), and a comma-separated scan list. One row per (session, line). Only rows matching `--line` are processed; each is loaded from a directory named `{Project}_{Session}` (e.g. `AGBT25B_386_03`) relative to the current directory.
- Per row, per feed: `sdf.getfs(...)` (frequency-switched, `ifnum=0, plnum=0`, channel-sliced, opacity-corrected, units Ta or Jy) is appended to `LineParams.final_sb`. If baseline order >= 0, `baseline()` fits a polynomial to each feed's time-average and subtracts it from the ScanBlock. `LineParams.write` concatenates feed ScanBlocks and writes SDFITS (default `pelican_{LINE}_Session{N}.sdfits`).
- Known quirks (verify before relying on them):   `--split` and `--nchan` are parsed but unused. Helper functions at the top (`plotbeams`, `ralongmapview`, `pointing_errors`, `specplot`) are for interactive/notebook use and rely on the `pelicancenter` global defined under `__main__`.
- Requires Python >= 3.10 (uses `str | u.Unit` annotations) and `dysh`.

## Batch driver: `doall`

Bash script. Per line (HCO+, HCN, CS) runs a pipeline, all lines concurrently: `pelican.py` per session in parallel -> (`-g`) one `gbtgridder` call over all that line's SDFITS, options from `grid.opts` + `--grid-args`, restfreq from `input.tab`, never `--scans` (scan numbers overlap across sessions) -> (`-m`) `domoment.py` (argparse wrapper around the `maskmoment` clone) on `pelican_{LINE}_cube.fits` only, options from `moment.opts` + `--moment-args`. Also `--dry-run`, `--verbose`. `gbtgridder/` and `maskmoment/` are git-ignored add-on clones.

## Observing-session housekeeping (Makefile + GBO scripts)

The Makefile is largely inherited from GBT-EDGE; many targets (`bench*`, `stats`, `all`, `NGC0001`, `gbtpipe`/`degas` installs) are unused leftovers. The ones that matter here (`make help` lists the `##` ones), steered by `SEQ` (session number, default from an optional `SEQ` file) and `PID=AGBT25B_386`:

- `make sdfits` / `make sdfits2`: list raw sessions / size of one in `/home/sdfits` (GBO only).
- `make rsync [SEQ=03]`: copy raw data and weather coefficients to Maryland.
- `make tsys [SEQ=03]`: runs `tsys.py` (GBTIDL-based, GBO only, slow), copies `.tsys` to `tsyslogs/` and `.summary` to `summarylogs/`; prints the `git add`/`commit` commands rather than running them.
- `make astrid [SEQ=03]`: runs `getastridlog` into `astridlogs/`.

`tsys.py` drives the IDL procedures in `pro/` (`vanecal*.pro`, `argus_*.pro`, `getatmos.pro`). `tsys_dysh.py` is a dysh-based alternative that computes VANE Tsys per scan and beam. `get_weather` is a bash loop that downloads the GBO weather forecast image periodically. See `README_workflow.md` for the observing workflow.

## Data files

`test.sdfits`, `input.tab` edits, and `.virtual_documents/` are local working files, not tracked in git (large SDFITS output especially; do not commit).
