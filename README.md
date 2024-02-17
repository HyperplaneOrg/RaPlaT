# RaPlaT

This code is a derivative of the Radio Planning Tool (RaPlaT) for the GRASS GIS System.
See the original COPYING file for additional information on copyrights and credits.
As of Feb 2024, a version of the original code can be found here:

[https://e6.ijs.si/RaPlaT/GRASS-RaPlaT_main_page.htm](https://e6.ijs.si/RaPlaT/GRASS-RaPlaT_main_page.htm)


## Building and Installing the Extensions Manually

Determine where this GRASS module source is located on the file system. Then, start grass and attempt to
build and install the extensions:

```shell
$ grass
GRASS : ~ > g.extension extension=r.raplat url=/Users/foobar/RaPlaT/python3/r.raplat
Compiling...
Installing...
Updating extensions metadata file...
Updating extension modules metadata file...
Installation of <r.raplat> successfully finished
GRASS : ~ >
```

Continue to Install the helper extensions:

```shell
GRASS : ~ > g.extension extension=r.MaxPower url=/Users/foobar/RaPlaT/raster/r.MaxPower
GRASS : ~ > g.extension extension=r.clutconvert url=/Users/foobar/RaPlaT/raster/r.clutconvert
GRASS : ~ > g.extension extension=r.cost231 url=/Users/foobar/RaPlaT/raster/r.cost231
GRASS : ~ > g.extension extension=r.fspl url=/Users/foobar/RaPlaT/raster/r.fsp
GRASS : ~ > g.extension extension=r.hata url=/Users/foobar/RaPlaT/raster/r.hata
GRASS : ~ > g.extension extension=r.hataDEM url=/Users/foobar/RaPlaT/raster/r.hataDEM
GRASS : ~ > g.extension extension=r.sector url=/Users/foobar/RaPlaT/raster/r.sector
GRASS : ~ > g.extension extension=r.waik url=/Users/foobar/RaPlaT/raster/r.waik
```

The example above works with GRASS 8.3.X on MacOS 13.X.X
One should be able to check the help:

```shell
GRASS : ~ > r.raplat --help
RaPlaT - raplat module (v30mar2023), RAdio PLanning Tool

Usage:
 r.raplat [-rpcx] csv_file=string [antmap_file=string] dem_map=string
   [clutter_map=string] [region=string] [rx_ant_height=value]
   [generate=string] [bandwidth=value] [rx_threshold=value] out_map=string
   [cellnum=value] [db_driver=string] [database=string] [out_table=string]
   [dbperf=value] [procnum=value] [freq_ovr=value] [radius_ovr=value]
   [model_ovr=string] [--overwrite] [--help] [--verbose] [--quiet]
   [--ui]

Flags:
  -r   Reuse results from existing intermediate model/sector files
  -p   (purge) Delete all unused sector radio coverage files
  -c   (check) Test run without actually performing radio coverage computation
  -x   Write all maps (model/sector/output) also as xyz files (in the working directory)

Parameters:
       csv_file   Radio cell/sector table in CSV format
    antmap_file   Antennas map file
                  default: $GISBASE/etc/radio_coverage/antenna_diagrams/antennamap.csv
        dem_map   DEM raster map for radio coverage simulation
    clutter_map   Clutter raster map (required for HataDEM model)
         region   Computation region (dem, current or region,raster,n,e,s,w,res (see g.region))
                  default: current
  rx_ant_height   Receiver antenna height [m]
                  default: 1.5
       generate   Selection of the generated output contents
                  options: rss-max,coverage,rss-sum,rss-maxix,lte-rssi,
                           lte-rsrp,lte-rsrq,lte-cinr,lte-maxspecteff,
                           lte-maxthrput,lte-interfere
                  default: rss-max
      bandwidth   Bandwidth [MHz]
                  default: 5
   rx_threshold   Minimum received power [dBm] for radio signal coverage
        out_map   Simulated radio coverage - raster (output)
                  default: out_raster
        cellnum   Number of succesive path loss values to be written in the table
                  default: 5
      db_driver   Database driver
                  options: none,dbf,mysql,pg,sqlite,csv
                  default: none
       database   Database name
                  default: $GISDBASE/$LOCATION_NAME/$MAPSET/dbf
      out_table   Simulated radio coverage - db table (output)
                  default: out_db
         dbperf   Database insert performance (rows/INSERT; 99: special fast mode via CSV)
                  options: 1-99
                  default: 20
        procnum   Number of parallel processes (-1: automatic, 0: non-parallel)
                  default: -1
       freq_ovr   Radio frequency override [MHz]
     radius_ovr   Radius override [km]
      model_ovr   Model override (with parameters)
```

Spot check the supporting modules/tools:

```shell
GRASS : ~ > r.fspl --help                                                                                                                       Mapset <PERMANENT> in <tmploc>
RaPlaT - Fspl module (v14aug2017)

Usage:
 r.fspl [-q] input_dem=name output=name [loss_exp=value]
   [loss_offset=value] coordinate=x,y [radius=value] [ant_height=value]
   [rx_ant_height=value] frequency=value [--overwrite] [--help]
   [--verbose] [--quiet] [--ui]

Flags:
  -q   Quiet

Parameters:
      input_dem   Name of input raster map
         output   Name for output raster map
       loss_exp   Exponent value (free-space-based fading model)
                  default: 2.0
    loss_offset   Offset value [dB] (free-space-based fading model)
                  default: 0.0
     coordinate   Base station coordinates
         radius   Computation radius [km]
                  default: 10
     ant_height   Transmitter antenna height [m]
                  default: 10
  rx_ant_height   Receiver antenna height [m]
                  default: 1.5
      frequency   Frequency [MHz]
GRASS : ~ >
```

## The RaPlaT Extension Documents

One can find the documents and manuals here:
[https://github.com/HyperplaneOrg/RaPlaT/tree/master/docs](https://github.com/HyperplaneOrg/RaPlaT/tree/master/docs)

