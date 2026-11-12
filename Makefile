#	 various target to aid in installing GBT-EDGE
#
SHELL = /bin/bash

# use standard wget or Peter's caching wgetc 
WGET = wget

# use python3 or anaconda3
PYTHON = anaconda3

# timing
TIME = /usr/bin/time

# parallel?
OMP = OMP_NUM_THREADS=1

# Greenbank locations (useful for rsync)

SDIR = /home/sdfits
WDIR = /users/rmaddale/Weather/ArchiveCoeffs
PID  = AGBT25B_386


# git directories we should have here and keep updated

GIT_DIRS = gbtpipe degas maskmoment gbtgridder

# URLs that we may need

URL0a = https://www.astro.umd.edu/~teuben/edge/data/AGBT21B_024_01.tar
URL0b = https://www.astro.umd.edu/~teuben/edge/data/GBTWeather.tar.gz
URL1  = https://github.com/GBTSpectroscopy/gbtpipe
URL1a = https://github.com/teuben/gbtpipe
URL2  = https://github.com/GBTSpectroscopy/degas
URL2a = https://github.com/teuben/degas
URL3  = https://github.com/astroumd/lmtoy
URL4  = https://github.com/tonywong94/maskmoment
URL4a = https://github.com/teuben/maskmoment
URL5  = https://github.com/tonywong94/edge_pydb
URL6  = https://github.com/richteague/bettermoments
URL7  = https://github.com/GreenBankObservatory/gbtgridder
URL7a = https://github.com/teuben/gbtgridder
URL7b = https://github.com/GreenBankObservatory/dysh # we need the latest version with the FS fix
URL8  = https://github.com/radio-astro-tools/spectral-cube/
URL9  = https://github.com/teuben/nemo

.PHONY:  help install build

## help:	 This Help
help : Makefile
	@sed -n 's/^##//p' $<

## git:		 make the needed git repos (GIT_DIRS)
git:  $(GIT_DIRS)
	@echo Last git: `date` >> git.log

## pull:	 update (git pull) the used git repos
pull:
	@echo -n "lmtoy: "; git pull
	-@for dir in $(GIT_DIRS); do\
	(echo -n "$$dir: " ;cd $$dir; git pull); done
	@echo Last pull: `date` >> git.log

## status:	 query status of the used git repos
status:
	@echo -n "lmtoy: "; git status -uno
	-@for dir in $(GIT_DIRS); do\
	(echo -n "$$dir: " ;cd $$dir; git status -uno); done

gbtpipe:
	# git clone $(URL1)
	git clone -b twarm2 $(URL1a)
	(cd $@; git remote add upstream $(URL1))

gbtpipe_pjt:
	git clone $(URL1a) gbtpipe_pjt
	(cd gbtpipe_pjt; git remote add upstream $(URL1))

degas:
	#git clone $(URL2)
	git clone -b tmb $(URL2a)
	(cd $@; git remote add upstream $(URL2))

degas_pjt:
	git clone $(URL2a) degas_pjt
	(cd degas_pjt; git remote add upstream $(URL2))

maskmoment:
	git clone $(URL4a)

dysh:
	git clone $(URL7b) dysh
	(cd dysh; pip install -e .)

gbtgridder_gbo:
	git clone $(URL7) gbtgridder_gbo

gbtgridder:
	# git clone -b release_3.0	$(URL7)
	git clone -b python3 $(URL7a)
	(cd gbtgridder; git remote add upstream $(URL7))  

spectral-cube:
	git clone $(URL8)

pel_env:
	python3 -m venv pel_env

install_anaconda:
	./install_anaconda3


# note gbtpipe needs to be installed before degas
install_gbtpipe:  gbtpipe pel_env
	(cd gbtpipe;\
	source ../pel_env/bin/activate;\
	pip3 install --upgrade pip;\
	pip3 install -e .)

install_degas:	degas pel_env
	(cd degas;\
	source ../pel_env/bin/activate;\
	pip3 install --upgrade pip;\
	pip3 install -e .)

install_maskmoment:  maskmoment pel_env
	(cd maskmoment;\
	source ../pel_env/bin/activate;\
	pip3 install --upgrade pip;\
	pip3 install -e .)

#  running at GBO, rawdata just points to /home/sdfits
#  offsite you will need to supply your own, YMMV
#  usually using the rsync target here, and running it from GBO
#  this also updates the weather database

## sdfits:	 show sessions in /home/sdfits for this project
sdfits:
	@ls -dl /home/sdfits/$(PID)_*


## sdfits2:  show filesize in selected SEQ
sdfits2:
	@ls -l /home/sdfits/$(PID)_$(SEQ)

## rawdata:  symlink to the SDFITS data (SDIR=)
rawdata:
	@if [ -d $(SDIR) ]; then \
	  ln -s $(SDIR) rawdata; \
	else \
	  echo "Not at GBO; Provide your own symlink/directory named 'rawdata'"; \
		  echo "or get a precomputed dataset via: make NGC0001"; \
	  echo "SDIR=$(SDIR)"; \
	fi

## weather:  symlink to the GBO weather database directory (WDIR=)
weather:
	@if [ -d $(WDIR) ]; then \
	  ln -s $(WDIR) weather; \
	else \
	  echo "Not at GBO; Provide your own symlink/directory named 'weather'"; \
	  echo "WDIR=$(WDIR)"; \
	fi

weather0:
	wget -q $(URL0b) -O - | tar zxf -
	ln -s GBTWeather weather


other:
	pip install radio_beam
	pip install spectral_cube

#  rsync:	only useful if you are at GBO and want to rsync to your offline location
#  if you want to fully reduce the data off-line
SEQ = 01
-include SEQ
REM = teuben@lma.astro.umd.edu:/n/lma1/teuben/
## rsync:	 rsync to REM=$REM and SEQ=$SEQ
rsync:
	@echo rsync to REM=$(REM) and SEQ=$(SEQ)
	du -sh $(SDIR)/$(PID)_$(SEQ)
	@echo weather
	-rsync -ahv --bwlimit=8000 $(WDIR)/Coeffs* $(REM)/GBTWeather
	@echo rawdata SEQ=$(SEQ)
	-rsync -ahv --bwlimit=8000 $(SDIR)/$(PID)_$(SEQ) $(REM)/GBTRawdata

#  this lengthy IDL based procedure computes the mean/rms/min/max for tsys for a given SEQ
## tsys:	 make tsys and summary files for SEQ=$SEQ (extra tsys.py options via TSYSARGS, e.g. TSYSARGS="--no-beams --vc3")
tsys:
	./tsys.py $(TSYSARGS) $(PID)_$(SEQ)
	mv pro/$(PID)_$(SEQ).tsys tsyslogs
	@echo git add			tsyslogs/$(PID)_$(SEQ).tsys
	@echo git commit -m new tsyslogs/$(PID)_$(SEQ).tsys
	mv pro/$(PID)_$(SEQ).summary summarylogs
	@echo git add			summarylogs/$(PID)_$(SEQ).summary
	@echo git commit -m new summarylogs/$(PID)_$(SEQ).summary

## astrid:	 make astridlogs for SEQ=$SEQ
astrid:
	(cd astridlogs; getastridlog $(PID)_$(SEQ))
	@echo git add			astridlogs/$(PID)_$(SEQ)_log.txt
	@echo git commit -m new astridlogs/$(PID)_$(SEQ)_log.txt

