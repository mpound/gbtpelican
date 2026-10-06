# Workflow

## FastX for VNC

https://ssh.gb.nrao.edu:3443

where you start astrid and cleo.


## Data view

During observing at GBT the raw data is written to /home/sdfits/AGBT25B_386_*

The command

      make sdfits

will show which ones exist, and their last modification date. The command

      make sdfits2

shows the size of the last (`SEQ`) session.

## Sync the data to maryland

During the observing you can rsync the data to Maryland, but you can also wait until the end
of course. Depending on data reduction scripts it may be useful to reduce on a much faster
Maryland machine.

It can only be done if you have write permission to /lma1/teuben/GBTRawdata/ (raw data) and
/lma1/teuben//GBTWeather (weather). The command

      make rsync [SEQ=03]

will do this.  Change the `SEQ` file if you don't want to repeat having to give the `SEQ=` argument.

## Final closeup

WHen all observing is done, the tsys, summary, and astridlogs should be created. It will remind
(via copy/paste) how to add them to this git repo so they can be pushed from GBO. No rsync here.

      make tsys [SEQ=03]

this command will take a while. It will create both the tsys and summary files.  The astrid logs with

      make astrid [SEQ=03]

will be much quicker.

