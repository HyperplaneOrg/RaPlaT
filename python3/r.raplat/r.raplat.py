#!/usr/bin/env python3

############################################################################
#
# MODULE:      r.raplat
# AUTHOR(S):   Igor Ozimek, Jozef Stefan Institute
#
# PURPOSE:     RAdio PLanning Tool - main script
#
# COPYRIGHT:   (C) 2015-2023 Jozef Stefan Institute
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with RaPlaT
#              for details.
#
#############################################################################
#
#
# Revision history r.radcov (OLD):
#
# 22-jan-10:
#  - during model/sector simulation display sector name (userLabel) and number of current/all simulated
#    models or sectors
#  - changed model and sector filenames, sector names now include (start with) the sector name (userLabel)
#    model filenames start with a leading underline
#
# 26-jan-10:
#  - options not used in the "frequency" parameter definition (limits allowed values to options value only)
#    default value remains 900 MHz
#  - db.GenerateTable and r.MaxPower not called if nothing to do (empty sector list)
#  - sector names as appearing in sector filenames are modified: any underline character is deleted
#    (this can possibly lead to duplicate sector names, but they are there only for human readibility)
#
# 29-jan-10:
#  - output table filename check according to the following rules:
#      supported table name characters are only: [A-Za-z][A-Za-z0-9_]*  (regular expression)
#      a table name must start with a character, not a number
#  - the above rule applies also to the output map file (just in case...)
#  - usrLabel string (sector name) is now also checked; new string type 'ul' is introduced for this purpose;
#    only teh following characters are allowed in usrLabel: 'A..Z', 'a..z', '0..9', '_', '-'
#  - modified allowed ranges (in cellTableDescrib):
#      electricalTiltAngle 0..10 (degrees)
#      power 0..50 (dBm - 1mW..100W)
#  - rules definition (cellTableDescrib) moved up, right after the module definition (grass.parser) part below
#
# 16-feb-10:
#  - radius value added to the model/sector filenames
#  - changed calling of (the new) r.sector - antennas folder parameter included (and default location changed)
#  - changed default input parameter values - CSV, DEM and clutter file names
#
# 24-feb-10:
#  - using new versions of db.CreateTable and r.MaxPower (supporting other databases besides dbf - e.g. MySQL)
#    additional input parameters for database driver and database (directory)
#
# 18-mar-10:
#  - modification for compatibility with GRASS 6.2 (originally the script was developed for GRASS 6.4)
#      appends path to grass.py module to python path (if not already there - as in the case of GRASS 6.4)
#      for GRASS 6.2: the modified grass.py module (the original is from GRASS 6.4) must be put in the same
#        directory as in GRASS 6.4 ($GISBASE/etc/python)
#  - changed default name and location of the antenna diagrams directory ($GISBASE/etc/antenna_diagrams)
#
# 9-avg-10  (a major new version):
#  - support for parallel (multicore) processing
#  - a new flag -c (check): test run only, without actually performing radio coverage computation
#  - command line parameter keywords changed according to GRASS standards
#    (lower case letters, underline as optional separatator)
#  - @MAPSET appended to GRASS input raster filenames for r.sector and r.MaxPower
#    (in the r.sector commands and in the sector file list for r.MaxPower)
#    this is to suppress possible warnings in case of duplicate model or sector filenames in other mapsets
#  - additional messages to report script execution progress
#  - times needed to compute models, sectors, and r.MaxPower are reported
#    (these are elapsed times, not the actual processing times - they depend on other running processes)
# The following changes require modifications to other radio coverage modules:
#  - use of the antennas map file to map antennas type and tilt to .MSI files
#    (instead of r.sector performing <antenna_type & el_tilt> to <.MSI_filename> translation by an internal algorithm)
#  - region management: computation region can be defined (set as temporal current region)
#    all coverage computations should be made within current region (instead of within input dem_map region),
#    requires changes to all coverage computation modules (models, r.sector, r.MaxPower)
#  - optional output data table creation (none if db_driver=_none_, modification of r.MaxPower is required) 
#
# 30-sep-10
#  - tested and corrected to actually work with the new versions of radio coverage computation modules
#  - additon of default_height and default_clutter parameters
#  - command line parameter value type checking
#    (grass parser does not check value types, e.g. for integer parameters)
#  - support (with automatic recognition) for the standard CSV format (delimiter ',') and
#    Microsoft Euro version (delimiter ';'=) (previuously only the MS Euro version was supported)
#  - region "refinement" (aligned to resolution, etc) for definiton options (previously only for current)
#  - dbdriver option _none_ changed to none
#
# 8-oct-10
#  - region "refinement" (alignment) algoritm changed (GRASS's built in g.region -a res=... doesn't work well)
#
# 15-oct-10
#  - cells within (extended) computation region but with transmitter location not defined in DEM are removed
#    from simulation, with warnings generated
#
# 22-oct-10
#  - corrected calling of r.sector: radius parameter missing - added
#  - "%" progress substrings postprocessing (removing) emproved (parallel processing only)
#  - added support for r.waik (Walfisch-Ikegami)
#  - updated/improved computation region management, 'dem' option added (default)
#
# 29-oct-10
#  - parameter P5 added to the models' (and sectors') filenames
#  - parameter P5 added to the temporary intermediate file for r.MaxPower
#
# 11-feb-11:
#  - modified for GRASS 6.4.0 (OS Ubuntu 10.04, Python 2.6.5)
#      changed module names
#      grass.region() now returns floatin point values (previously strings)
#
# 14-mar-11:
#  - modified for the new version of r.waik - clutter file not used any more
#
# 31-mar-11:
#  - modified for the newest version of r.waik - clutter file used again
#  - default_height & default_clutter not used any more
#
#
# 20-dec-11: (A major revison - modified for the new versions of RaPlaT C-language moduls (models, MaxPower))
#  - modifications for the r.MaxPower module (integrated db.GenerateTable funcionality, new parameter dbperf)
#  - some parameter modifications for existing models (e.g. waik does not use clutter file any more, ...)
#  - support for additional models fspl and itm
#  - number of model parameters (input CSV file) enlarged from 5 to 11 (as required by the new version of waik model)
#      (processing generalized for possible future changes)
#  - parallel dispatch loop sleep time reduced to 10ms (from 0,2s, i.e. 20x)
#      estimated CPU overhead (idle dispatch loop) is <=0,2% of one CPU core
#  - enlarged position limits (in cellTableDescrib) to cover the entire Slovenia map
#  - userLabel replaced by cellName and antennaID (to support the new csv format
#      (written as cellName-antennaID in filenames)
#  - modified r.MaxPower: added ID field in output datatable (descriptive string fields retained for now)
#  - new override parameters: model_ovr and radius_ovr
#      (globally override radius and model w/ parameter values defined in radio cell/sector table for each cell/sector)
#  - grass.info() seems to fail for large messages (nothing printed) (commented aout and replaced by print in a few places)
#  - default region changed to 'current' (was 'dem')
#  - antennaType can now contain characters (in any order): alpha, num, ' ', '-', '/'
#      (previously it was limited to numbers only)
#  - -r now means the opposite(!): 'reuse' instead of 'recalculate' (default behaivour now is to recalculate all models/sectors)
#  - resolution added to the region report (previously only coordinates were printed)
#  - floating numbers allowed for region parameters values (previously integers only)
#  - region refinement (rounding of coordinates values according to the resolution) is now based on 'pixel centers'
#      (previously on 'pixel edges', i.e. region boundaries)
#      the improved algorithm works correctly also for negative coordinates values (the previous one did not)
#  - various minor changes and improvments
# 
# 25-jan-12:
#  - some minor bug corrections (including modifications of cellTableDescrib)
#
# 17-feb-12:
#  - maximum transmission power changed to 100 dBm (was 50 dBm)
#  - changed format of the intermediate (temporary) radio cell/sector table file for r.MaxPower: the first line (telling
#    the number of the following lines) is removed (r.MaxPower changed, does not need this any more)
# Possible future changes:
# - Individual transmission frequencies would be given for each cell in the input radio sector CSV file (instead globally
#   with the parameter 'frequnecy'
# - Antenna mapping: would not suppose .MSI / .msi file extensions, instead extesions would be explicitly given as part of
#   the filenames in the antenna mapping table
#
#
#  --- END OF DEVELOPMENT FOR MOBITEL/TELEKOM ---
#  (from now on they use their own modified versions of r.radcov)
#
#
# MODIFICATIONS FOR THE CREW PROJECT
# 18-oct-12:
#  - radius changed from integer to floating point (allowing for radius smaller then 1 km, with 1 m resolution)
#    (this was actually an error - it should already be floating point)
#  - radius_ovr is also changed to float
#  - power range is now -30..+100 dBm (was 0..+100dBm)
#  - additional parameter rx_threshold (optional, backward compatibility is manintained)
#    (this parameter is used by r.MaxPower)
#  - flag -1 added: Rx (dBm) values in output map replaced by 1.0 when above rx_threshold
#
# 13-nov-12:
#  - calling r.MaxPower_crew instead r.MaxPower (canceled 20-sep-13)
#
# 16-nov-12:
#  - added flag -x : all generated output maps (models/sectors/output) are output also as xyz files in the working directory
#
# 10-dec-12:
#  - purge list is displayed using print instead of grass_info. (grass_info() cannot handle very large strings, script crashes)
#
# 20-sep-13:
#  - merging with the IJS-E6 internal version (developed after 17-feb-12 until 20-sep-13):
#     - 14-feb-13:
#        - support added for r.urban model
#        - added parameter buildings_map (required for r.urban)
#        - removed the default value for the clutter_map parameter
#     - 01-avg-13:
#        - updated support for the new version of r.urban - additional parameter 'screen_separation'
#  - calling r.MaxPower again, instead of r.MaxPower_crew (undoing modification from 13-nov-12)
#  - ericsson changed to hataDEM
#
# 24-sep-13:
# - some clean-up of code
#
# 24-jul-14:
# - support for r.MaxPower quasi database type 'csv' added (this generates a CSV-format file instead of a database table)
# - added support for the receive antenna height command line parameter - currently only for hata, hataDEM and fspl propagation models
# - added support for the exponent and offset parameters for fspl (must be specified in the input csv table, cellTableDescrib was updated)
#
# 25-aug-14:
# - the parameter 'frequency' type changed to double (was integer)
# - all numeric command-line parameters (integer and double) are now checked for correctness (grass Python parser does not perform this)
# - added r.MaxPower parameter 'generate' (= rss-max, coverage, rss-maxix, rss-sum, lte-*...)
#   'rss-max' (default) generates the old output map, 'coverage' generates the old '-1' (flag) map
# - flag '-1' removed (superseded by 'generate=coverage')
# - added r.MaxPower parameter 'bandwidth' (actually required only for lte-*..., but not tested here)
# 
#
#------------------------------------------------------------------------------
# 
#
# Revision history r.raplat:
#
# 29-jun-15:
#  - based on r.radcov v25-aug-14 (currently its last version, no new version is planned)
#  - max antenna height changed to 20000m above the terrain (to allow for high altitude platforms)
#  - freq column added to the csv table
#  - column names and order of the csv table changed for better coherence
#  - parameter frequency changed to freq_ovr (since it is now actually frequency override)
#  - dropped compatibility with GRASS 6.2
#
# 16-aug-16
#  - support for the propagation model ITUR1546-4
#
# 16-aug-2017
# - modified for GRASS GIS 7 (was GRASS GIS 6)
#     option name in called RaPlaT modules changed: inputDEM -> input_dem; A0..A3 -> a0..a3, PHI_Street -> phi_street
#     changed parameters name in GRASS's r.what command
#     changed GRASS's command name g.mlist > g.list, an additional required parameter added
#
# 22-aug-2017
# - modified for GRASS GIS 7 (was GRASS GIS 6)
#     changed parameter name in GRASS's g.region command: rast -> raster
# - itm and urban model removed (code commented out)
# 
# 23-aug-2017
# - changed ITUR1546-4 option name rec_height -> rx_ant_height
#
# 24-aug-2017
# - added support for sqlite database (must be supported in the r.MaxPower module)
# - disabled support for the r.ITUR1546-4 module
# - additional correction for g.region command: rast -> raster
# 
# 10-dec-2018
#  - raster maps' color table changed to "rainbow" (which was default in GRASS 6, but not in GRASS 7)
#
# 6-jul-2020
# - removed #%Module "label" and modified "description" accordingly, for compatibility with GRASS 7.4
#
# 30.mar-2023
# - converted to Python3 (was Python2)
#
# Todo: reorganize python code (imports at the beginning, main, functions,,,)
#


#%Module
#%  description: RaPlaT - raplat module (v30mar2023), RAdio PLanning Tool
#%end

#%option
#%  key: csv_file
#%  type: string
#%  gisprompt: old_file,file,input
#%  description: Radio cell/sector table in CSV format
##%  answer: sector_table.csv
#%  required: yes
#%end

#%option
#%  key: antmap_file
#%  type: string
#%  gisprompt: old_file,file,input
#%  description: Antennas map file
#%  answer: $GISBASE/etc/radio_coverage/antenna_diagrams/antennamap.csv
#%  required: no
#%end

#%option
#%  key: dem_map
#%  type: string
#%  gisprompt: old,cell,raster
#%  description: DEM raster map for radio coverage simulation
##%  answer: dem_map@PERMANENT
#%  required: yes
#%end

#%option
#%  key: clutter_map
#%  type: string
#%  gisprompt: old,cell,raster
##%  description: Clutter raster map (required for HataDEM and Urban models)
#%  description: Clutter raster map (required for HataDEM model)
##%  answer: clutter_map@PERMANENT
#%  required: no
#%end

##%option
##%  key: buildings_map
##%  type: string
##%  gisprompt: old,cell,raster
##%  description: Buildings raster map (required for Urban model)
###%  answer: buildings_map@PERMANENT
##%  required: no
##%end

#%option
#%  key: region
#%  type: string
#%  description: Computation region (dem, current or region,raster,n,e,s,w,res (see g.region))
#%  answer: current
#%  required: no
#%end

#%option
#%  key: rx_ant_height
#%  type: double
#%  description: Receiver antenna height [m]
#%  answer: 1.5
#%  required: no
#%end

#%option
#%  key: generate
#%  type: string
#%  description: Selection of the generated output contents
#%  options: rss-max,coverage,rss-sum,rss-maxix,lte-rssi,lte-rsrp,lte-rsrq,lte-cinr,lte-maxspecteff,lte-maxthrput,lte-interfere
#%  answer: rss-max
#%  required: no
#%end

#%option
#%  key: bandwidth
#%  type: double
#%  description: Bandwidth [MHz]
#%  answer: 5
#%  required: no
#%end

#%option
#%  key: rx_threshold
#%  type: double
#%  description: Minimum received power [dBm] for radio signal coverage
#%  required: no
#%end

#%option
#%  key: out_map
#%  type: string
#%  gisprompt: new,cell,raster
#%  description: Simulated radio coverage - raster (output)
#%  answer: out_raster
#%  required: yes
#%end

#%option
#%  key: cellnum
#%  type: integer
#%  description: Number of succesive path loss values to be written in the table
#%  answer: 5
#%  required: no
#%end

#%option
#%  key: db_driver
#%  type: string
#%  description: Database driver
#%  options: none,dbf,mysql,pg,sqlite,csv
##%  options: none,dbf,mysql,odbc,ogr,pg,sqlite
#%  answer: none
#%  required: no
#%end

#%option
#%  key: database
#%  type: string
#%  description: Database name
#%  answer: $GISDBASE/$LOCATION_NAME/$MAPSET/dbf
#%  required: no
#%end

#%option
#%  key: out_table
#%  type: string
#%  gisprompt: new_dbtable,dbtable,dbtable
#%  description: Simulated radio coverage - db table (output)
#%  required: yes
#%  answer:out_db
#%  required: no
#%end

#%option
#%  key: dbperf
#%  type: integer
#%  description: Database insert performance (rows/INSERT; 99: special fast mode via CSV)
#%  options: 1-99
#%  answer: 20
#%  required: no
#%end

#%option
#%  key: procnum
#%  type: integer
#%  description: Number of parallel processes (-1: automatic, 0: non-parallel)
#%  answer: -1
#%  required: no
#%end

#%option
#%  key: freq_ovr
#%  type: double
#%  description: Radio frequency override [MHz]
#%  required: no
#%end

#%option
#%  key: radius_ovr
#%  type: double
#%  description: Radius override [km]
#%  required: no
#%end

#%option
#%  key: model_ovr
#%  type: string
#%  description: Model override (with parameters)
#%  required: no
#%end

#%flag
#%  key: reuse
#%  description: Reuse results from existing intermediate model/sector files
#%end

#%flag
#%  key: purge
#%  description: (purge) Delete all unused sector radio coverage files
#%end

#%flag
#%  key: check
#%  description: (check) Test run without actually performing radio coverage computation
#%end

#%flag
#%  key: x
#%  description: Write all maps (model/sector/output) also as xyz files (in the working directory) 
#%end


# ---- RULES DEFINITION FOR THE INPUT (CSV) TABLE ----

cellTableDescrib = [
                    ['cellName', 'name'],
                    ['antID', 'id', [1, 999999]],
                    ['antType', 'antype'],
                    ['antEast', 'i', [0, 0]],  # any region (Slovenia: [370000, 630000])
                    ['antNorth', 'i', [0, 0]],  # any region (Slovenia: [25000, 200000])
                    ['antHeightAG', 'f', [0., 20000.]], # heigt above ground
                    ['antDirection', 'i', [0, 360]],
                    ['antElecTilt', 'i', [0, 10]],
                    ['antMechTilt', 'i', [-90, +90]],
                    ['freq', 'f', [1., 10000.]],
                    ['power', 'f', [-30., 100.]],
                    ['radius', 'f', [0., 100.]],
                    ['model', 's', ['hata'],['cost231'],['hataDEM'],['waik'],['fspl']],
                    ['P1', 's', ['urban', 'suburban', 'open'],\
                           's', ['metropolitan', 'medium_cities'],\
                           'f', [0., 0.],\
                           'f', [20., 60.],\
                           'f', [0., 9.]],
#                           hata cost231 hataDEM        waik                                    fspl
                    [ 'P2', '-', '-',    'f', [0., 0.], 'i', [30,  70],                         'f', [-100., +100.]], 
                    [ 'P3', '-', '-',    'f', [0., 0.], 'i', [ 5,  35],                         '-'],
                    [ 'P4', '-', '-',    'f', [0., 0.], 'i', [ 3,  15],                         '-'],
                    [ 'P5', '-', '-',    '-',           'i', [ 3,  25],                         '-'],
                    [ 'P6', '-', '-',    '-',           'i', [10,  30],                         '-'],
                    [ 'P7', '-', '-',    '-',           'i', [10,  25],                         '-'],
                    [ 'P8', '-', '-',    '-',           'i', [20,  50],                         '-'],
                    [ 'P9', '-', '-',    '-',           'i', [ 0, 300],                         '-'],
                    ['P10', '-', '-',    '-',           'i', [ 0, 180],                         '-'],
                    ['P11', '-', '-',    '-',           's', ['metropolitan', 'medium_cities'], '-']
                   ]
#
#                    ['model', 's', ['hata'],['cost231'],['hataDEM'],['waik'],['fspl'],['itm'],['urban'],['ITUR1546-4']],
#                    ['P1', 's', ['urban', 'suburban', 'open'],\
#                           's', ['metropolitan', 'medium_cities'],\
#                           'f', [0., 0.],\
#                           'f', [20., 60.],\
#                           'f', [0., 9.],\
#                           'f', [0., 0.],\
#                           's', ['half-screen', 'microcell', 'half-screen+microcell'],\
#                           's', ['urban', 'suburban', 'rural', 'none']],
##                           hata cost231 hataDEM        waik                                    fspl                 itm               urban            ITUR1546-4
#                    [ 'P2', '-', '-',    'f', [0., 0.], 'i', [30,  70],                         'f', [-100., +100.], 'f', [0., 0.],    'i', [1, 10000], 'f', [0., 100.]], 
#                    [ 'P3', '-', '-',    'f', [0., 0.], 'i', [ 5,  35],                         '-',                 'f', [0., 0.],    '-',             '-'],
#                    [ 'P4', '-', '-',    'f', [0., 0.], 'i', [ 3,  15],                         '-',                 'i', [0, 0],      '-',             '-'],
#                    [ 'P5', '-', '-',    '-',           'i', [ 3,  25],                         '-',                 'i', [0, 1],      '-',             '-'],
#                    [ 'P6', '-', '-',    '-',           'i', [10,  30],                         '-',                 'f', [0.1, 0.99], '-',             '-'],
#                    [ 'P7', '-', '-',    '-',           'i', [10,  25],                         '-',                 'f', [0.1, 0.99], '-',             '-'],
#                    [ 'P8', '-', '-',    '-',           'i', [20,  50],                         '-',                 '-',              '-',             '-'],
#                    [ 'P9', '-', '-',    '-',           'i', [ 0, 300],                         '-',                 '-',              '-',             '-'],
#                    ['P10', '-', '-',    '-',           'i', [ 0, 180],                         '-',                 '-',              '-',             '-'],
#                    ['P11', '-', '-',    '-',           's', ['metropolitan', 'medium_cities'], '-',                 '-',              '-',             '-']
#                   ]


import os, sys
if "GISBASE" not in os.environ:
    print("You must be in GRASS GIS to run this program.")
    sys.exit(1)


#import grass
from grass.script import core as grass
options, flags = grass.parser()

mapset = grass.gisenv()['MAPSET']
#workDir = grass.gisenv()['GISDBASE'] + '/' + grass.gisenv()['LOCATION_NAME'] + '/' + mapset


# check output table and map filenames
import re

outFilename = options['out_map']
if re.match(r'[A-Za-z][A-Za-z0-9_]*\Z', outFilename) == None:
    grass.fatal("Wrong output map filename '" + outFilename + "'\n" +
                "allowed chars: 'A-Z', 'a-z', '0-9', and '_' (must start with a char)")

dbFilename = options['out_table']
if re.match(r'[A-Za-z][A-Za-z0-9_]*\Z', dbFilename) == None:
    grass.fatal("Wrong output table filename '" + dbFilename + "'\n" +
                "allowed chars: 'A-Z', 'a-z', '0-9', and '_' (must start with a char)")


# check numeric arguments values (grass parser does not perform this)
intArgList = ['cellnum', 'dbperf', 'procnum']
optIntArgList = []
floatArgList = ['rx_ant_height', 'bandwidth']
optFloatArgList = ['freq_ovr', 'radius_ovr', 'rx_threshold']
for optKey in ( intArgList + optIntArgList + floatArgList + optFloatArgList):
    optStr = options[optKey]
    if optKey in ( optIntArgList + optFloatArgList):
        if optStr == '':
            continue
    if optKey in ( intArgList + optIntArgList):
        try:
            int(optStr)
        except:
            grass.fatal('Non-integer value for integer parameter: ' + optKey + '=' + optStr )
    else:
        try:
            float(optStr)
        except:
            grass.fatal('Non-numeric value for float parameter: ' + optKey + '=' + optStr )


# ---- SET THE NUMBER OF PARALLEL PROCESSES
# ---- used for parallel computations of model and sectors

# detectCPU is from:
#   http://www.artima.com/weblogs/viewpost.jsp?thread=230001
#   Computing Thoughts
#   Concurrency with Python, Twisted, and Flex
#   by Bruce Eckel
#   May 3, 2008

def detectCPUs():
    """
    Detects the number of CPUs on a system. Cribbed from pp.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if "NUMBER_OF_PROCESSORS" in os.environ:
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
            if ncpus > 0:
                return ncpus
    return 1 # Default

#----------------

parNum = int(options['procnum'])
if parNum < 0:
    parNum = detectCPUs()
    grass.info('\nNumber of detected processors = ' + str(parNum))


# ---- READ IN THE ANTENNAS-MAPPING CSV TABLE ----

grass.info('\nREADING AND CHECKING THE ANTENNAS-MAPPING CSV TABLE ...') 

import csv

# get the name of the antenna map file from the parameter (or default)
amPathname = options['antmap_file']
# evaluate GISBASE enviroment variable
gisBase = os.getenv('GISBASE')
amPathname = amPathname.replace('$GISBASE',gisBase)
# create full absolute path for inPathname
amPathname = os.path.join(gisBase + '/etc/radio_coverage/antenna_diagrams', amPathname)

try:
    csvFile = open(amPathname, 'r')
except IOError as xxx_todo_changeme:
    (errno, strerror) = xxx_todo_changeme.args
    grass.fatal("Cannot open Antennas-Mapping input file '" + amPathname + "'")
    
#standard type of CSV required (delimiter ',')
csvReader = csv.reader( csvFile)

antMapTable=[]
for row in csvReader:
    if row == []:  # skip empty lines
        continue
    if row[0][0] != '#':  # append line if not commented out
        antMapTable.append(row)    
csvFile.close()


# ---- CHECK AND REMOVE THE HEADER (the first line) ----

refAntMapTableHeader = ['antennaType', 'frequency', 'frequencyLower', 'frequencyUpper', 'EDT', 'MSIfilename', 'technology']
antMapTableHeader = antMapTable[0]
if antMapTableHeader != refAntMapTableHeader:
    _msgStr = 'Error in Antennas-Mapping Table header:\n'
    _msgStr += str(antMapTableHeader)
    _msgStr += '\nShould be:\n'
    _msgStr += str(refAntMapTableHeader)
    grass.fatal(_msgStr)
del antMapTable[0]


# ---- CHECK COLUMN VALUES, BUILD .MSI FILES DICTIONARY (short filename -> full pathname) ----
# Column value rules:
# antennaType: [alpha digit - / .] (in any order)
# frequency (_Lower, _Upper): integer/floating, frequencyLower < frequency < frequencyUpper
# EDT: integer/floating, 0..90
# MSIFIleName: extension should be .MSI, file located in a subdriectory, no duplicates allowed (same filenames in different directories)
#              given (bare) filenames are replaced with full paths (with extension .MSI)
# technology: currently [GSM900 GSM1800 UMTS2100], but not checked (ignored)
#
# table row duplicates (same antennaType & frequency band & EDR w/ possibly different MSIFileName) are not checked
#   (first suitable antenna diagram is used)

# make a list of all files in the antenna_diagrams subdirectories, check for duplicates, create dictionary (name -> fullpath)
amDir = os.path.dirname(amPathname)  # the antenna_daigrams directory (location of the antennamap file)
amFileList = []
errf = False
for rootName, dirList, fileList in os.walk(amDir):
    for filename in fileList:
        [fnBase, fnExt] = os.path.splitext(filename)
        if fnExt.upper() == '.MSI':
            if fnBase in [amFileListRow[0] for amFileListRow in amFileList]:
                grass.error("Duplicate .MSI file found in the directory tree: '" + filename + "'")
                errf = True
            else:
                amFileList.append([fnBase, os.path.join(rootName, filename)])
amDict = dict(amFileList)  # create filenames dictionary
if errf:
    sys.exit(1)
grass.info('Number of .MSI files found: ' + str(len(amFileList)))

# Check column values
import string
errmsg = ''
antTypeChars = string.ascii_letters + string.digits + ' ' + '-' + '/' + '.'
for antMapTableRow in antMapTable:
    [antType, antFreq, antFreqLow, antFreqHigh, antTilt, MSIfname, tech] = antMapTableRow
    fAntFreq     = float(antFreq)
    fAntFreqLow  = float(antFreqLow)
    fAntFreqHigh = float(antFreqHigh)
    fAntTilt     = float(antTilt)
    for atchar in antType:
        if not atchar in antTypeChars:
            errmsg += 'Illegal caharacters in Antenna type, allowed: ' + antTypeChars + '\n'
            break
    if fAntFreqLow > fAntFreqHigh:
        errmsg += 'Frequency error: Lower > Upper\n'
    if fAntFreq < fAntFreqLow:
        errmsg += 'Frequency error: < Lower\n'
    if fAntFreq > fAntFreqHigh:
        errmsg += 'Frequency error: > Upper\n'
    if fAntTilt < 0 or fAntTilt > 90:
        errmsg += 'EDT of of range, allowed: 0..90\n'
    if not MSIfname in list(amDict.keys()):
        errmsg += "MSI file not found: '" + MSIfname + "'\n"
    if errmsg != '':
        break
if errmsg != '':
    grass.fatal(errmsg + ' in row:\n' + str( antMapTableRow))


def findAntennas( aMapTable, rType, rFreqNum, rTiltNum):
    """ Find suitable antennas (returns a list) """
    antennas = []
    for [aType, aFreq, aFreqLow, aFreqHigh, aTilt, aMSIfname, aTech] in aMapTable:
        if (rType == aType and
            rFreqNum >= float(aFreqLow) and rFreqNum <= float(aFreqHigh) and
            rTiltNum == float(aTilt)):
            antennas.append( (aMSIfname, abs(rFreqNum - float(aFreq))) )
    if len(antennas) == 0:
        return []
    antennas.sort( key=lambda row: row[1])  # sort by frequency difference (nominal freq.)
    return [row[0] for row in antennas]


# ---- READ IN THE RADIO SECTOR CSV TABLE ----

grass.info('\nREADING AND CHECKING THE RADIO SECTOR CSV TABLE ...') 

#create reference header from cellTableDescrib
#refCellTableHeader = ['cellName', 'antID', 'antType', ...
refCellTableHeader = []
for clmnlist in cellTableDescrib:
    refCellTableHeader.append(clmnlist[0])

#open the radio sector CSV file
try:
    csvFile = open(options['csv_file'], 'r')
except IOError:
    grass.fatal('Cannot open Radio cell/sector input file')
    
#determine the type of CSV - standard (delimiter ',') or MS Euro (delimiter ';')
hline1 = csvFile.readline()
hline1 = re.sub(r'".*?"', '', hline1)  #remove quoted text
hStd = hline1.count(',') > 0
hMSE = hline1.count(';') > 0
if not hStd and not hMSE:
    grass.fatal('CSV file first line (header) error: no separator (commas or semicolons) found')
if hStd and hMSE:
    grass.fatal('CSV file first line (header) error: separator ambiguity - commas and semicolons found')
if hStd:
    grass.info('Standard CSV file format detected')
    csvReader = csv.reader( csvFile)
else:
    grass.info('MS Euro CSV file format detected')
##    csvReader = csv.reader( csvFile, delimiter=';', quoting=csv.QUOTE_NONE)
    csvReader = csv.reader( csvFile, delimiter=';')
csvFile.seek(0)

cellTable=[]
for row in csvReader:
    if row == []:  # skip empty lines
        continue
    if row[0][0] != '#':  # append line if not commented out
        for ix in range( len(row), len(refCellTableHeader)):
            row.append('')   #add dummy items to the row until reference header length (number of columns) is reached
        cellTable.append(row)    

csvFile.close()


# ---- CHECK AND REMOVE THE HEADER (the first line) ----

# extract and check the header (and delete it from the table)
cellTableHeader = cellTable[0]
if cellTableHeader != refCellTableHeader:
    _msgStr = 'Error in Cell Table header:\n'
    _msgStr += str(cellTableHeader)
    _msgStr += '\nShould be:\n'
    _msgStr += str(refCellTableHeader)
    grass.fatal(_msgStr)
del cellTable[0]


# ---- OVERRIDE PROCESSING ----
# freq, radius and model w/ parameters given in radio sector csv table can be overriden

overrideList = []
freqOvrStr = options['freq_ovr']
radiusOvrStr = options['radius_ovr']
modelOvrStr = options['model_ovr']

if freqOvrStr != '':
    overrideList.append( ('freq', freqOvrStr) )

if radiusOvrStr != '':
    overrideList.append( ('radius', radiusOvrStr) )

if modelOvrStr != '':
    modelOvrList = modelOvrStr.split(',')
    overrideList.append( ('model', modelOvrList[0])) 
    ix = -1
    for ix,paramValueStr in enumerate(modelOvrList[1:]):
        overrideList.append( ('P' + str(ix+1), paramValueStr))
    ix += 2
    while 'P' + str(ix) in refCellTableHeader:  # clear the rest of Pn's
        overrideList.append( ('P' + str(ix), ''))
        ix += 1

for ovrItem in overrideList:
    if ovrItem[0] in refCellTableHeader:
        clmnix = refCellTableHeader.index(ovrItem[0])
        for rowix,row in enumerate(cellTable):  # step through cellTable lines
                cellTable[rowix][clmnix]=ovrItem[1]  # replace column
    else:
        print('??? INTERNAL ERROR (overrideList) - override item not in refCellTableHeader')
        sys.exit(1)


# ---- CHECK 'cellName' AND 'antID' COLUMNS FOR DUPLICATES (warning and error respectivelly)

for tcix, clmnDescr in enumerate(cellTableDescrib):  # tcix: cellTable column index
    if clmnDescr[1] in ['id','name']:
        dupCheckList = []
        for cellTableLine in cellTable:
            cellTableField = cellTableLine[tcix]
            if cellTableField in dupCheckList:
                _msgStr = " Duplicated '" + clmnDescr[0] + "': " + cellTableField
                if clmnDescr[1] == 'id':
                    grass.fatal(_msgStr)
                else:
                    pass  # no warnings needed for duplicate names (it is normal)
##                    grass.warning(_msgStr)
            else:
                dupCheckList.append(cellTableField)
            

# ---- CHECK PARAMETERS (each table row) ----

# checking the input (CSV) table according to the rules given by cellTableDescrib
# a bit tricky and hard-to-explain procedural part of code follows :(
err = False  # becomes True if one or more parameter errors are found in a a table row
for cellTableLine in cellTable:
    vix = 1  # variant index (=1 unless multi-variant column definitions are used)
    for tcix, clmnDescr in enumerate(cellTableDescrib):  # tcix: cellTable column index
        # find the (chosen variant of the) parameter description list for the current parameter (table column)
        ibegin = 1  # this will be the beginning of the (chosen variant of the) parameter description list
        iend = len(clmnDescr)  # this will be the end of the (chosen variant of the) parameter description list
        n = 0  # a local variant counter (for a column)
        for i in range(1, iend):
            if isinstance(clmnDescr[i], str):
                n += 1
                if n == vix:
                    ibegin = i
                    break
        for i in range( ibegin + 1, iend):
            if not isinstance( clmnDescr[i], list):
                iend = i
                break

        # extract the (chosen variant of the) parameter description list for the current parameter
        clmnDescSublist = clmnDescr[ibegin:iend]


        # process (check correctness of) parameters according to cellTableDescrib
        if clmnDescSublist == []:
            print('??? INTERNAL ERROR (cellTableDescrib) - empty (or too short) column description')
            sys.exit(1)

        # unchecked type
        elif clmnDescSublist[0] == '-':
            if len(clmnDescSublist) != 1:
                print("??? INTERNAL ERROR (cellTableDescrib) - type '-' should have no additional parameters")
                sys.exit(1)

        # numeric types (integer, float) - check against min and max allowed values
        elif clmnDescSublist[0] in ['id', 'i', 'f']:
            if len(clmnDescSublist) != 2:
                print("??? INTERNAL ERROR (cellTableDescrib) - type 'i', 'id' or 'f' should have one additonal parameter (list)")
                sys.exit(1)
            nerr = (len(clmnDescSublist[1]) != 2)
            if not nerr:
                nummin = clmnDescSublist[1][0]
                nummax = clmnDescSublist[1][1]
            if clmnDescSublist[0] == 'f':
                nerr = nerr or not isinstance(nummin, float) or not isinstance(nummax, float)
                if nerr:
                    print('??? INTERNAL ERROR (cellTableDescrib) - invalid floating number range definition ( = ', clmnDescSublist[1], ')')
                    sys.exit(1)
            else:
                nerr = nerr or not isinstance(nummin, int) or not isinstance(nummax, int)
                if nerr:
                    print('??? INTERNAL ERROR (cellTableDescrib) - invalid integer number range definition ( = ', clmnDescSublist[1], ')')
                    sys.exit(1)

            # now process the value from the CSV table
            # but first change any (decimal) comma do (decimal) point
            #   (and store it back also to the original table)
            cellTableLine[tcix] = cellTableLine[tcix].replace(',', '.')
            _msgstr = 'Error in input table:\n' + ', '.join(cellTableLine) + '\n'
            try:
                if clmnDescSublist[0] == 'f':
                    num = float( cellTableLine[tcix])
                else:
                    num = int( cellTableLine[tcix])
            except ValueError:
                _msgstr += 'Column ' + str(tcix+1) + ' (=' + str(cellTableLine[tcix]) + ') should be '
                if clmnDescSublist[0] == 'f':
                    _msgstr += 'a floating point number'
                else:
                    _msgstr += 'an integer number'
                grass.error(_msgstr)
                err = True
            except IndexError:
                _msgstr +=  'Line too short - column ' + str(tcix+1) + ' missing (number expected)'
                grass.error(_msgstr)
                err = True
            if nummin < nummax:
                if (num < nummin) or (num > nummax):
                    _msgstr += ('Number in column ' + str(tcix+1) + ' (=' + str(cellTableLine[tcix]) + ') is out of range ' +
                                str(nummin) + '..' + str(nummax))
                    grass.error(_msgstr)
                    err = True
            
        # string type - check against allowed list of string - possibly more than one list (variants)
        elif clmnDescSublist[0] == 's':
            if len(clmnDescSublist) < 2:
                print("??? INTERNAL ERROR - type 's' should have one or more additonal parameter(s) (list(s))")
                sys.exit(1)
            ixstr = 0
            for strlst in clmnDescSublist[1:]:
                ixstr += 1
                if not isinstance(strlst, list):
                    print("??? INTERNAL ERROR (cellTableDescrib) - 's'  should be followed by lists(s) of strings")
                    sys.exit(1)
                # now process the value from the CSV table
                _msgstr = 'Error in input table:\n' + ', '.join(cellTableLine) + '\n'
                try:
                    if cellTableLine[tcix] in strlst:
                        vix *= ixstr  # matching string found - OK
                        break
                except IndexError:
                    _msgstr += 'Line too short - column ' + str(tcix+1) + ' missing (string expected)'
                    grass.error(_msgstr)
                    err = True
                    break
            else:
                _msgstr = 'Error in input table:\n' + ', '.join(cellTableLine) + '\n'
                _msgstr += 'String in column ' + str(tcix+1) + ' (=' + str(cellTableLine[tcix]) + ') is not allowed. Allowed: '
                fullstrlst = []
                for strlst in clmnDescSublist[1:]:
                    fullstrlst.extend(strlst)
                _msgstr += ', '.join(fullstrlst)
                grass.error(_msgstr)
                err = True

        # "name" string type - the following characters are allowed:
        #     'A..Z', 'a..z', '0..9', '_', '-'
        elif clmnDescSublist[0] == 'name':
            if len(clmnDescSublist) != 1:
                print("??? INTERNAL ERROR cellTableDescrib - type 'name' should have no additonal parameter (list)")
                sys.exit(1)
            if re.search(r'[^A-Za-z0-9_-]', cellTableLine[tcix]):
                _msgstr = 'Error in input table:\n' + ', '.join(cellTableLine) + '\n'
                _msgstr += 'Wrong string in column ' + str(tcix+1) + ' (=' + str(cellTableLine[tcix]) + ').'
                _msgstr += " Allowed chars: 'A..Z', 'a..z', '0..9', '_', '-'"
                grass.error(_msgstr)
                err = True

        # "antype" string type - the following characters are allowed:
        #     'A..Z', 'a..z', '0..9', '_', '-'
        elif clmnDescSublist[0] == 'antype':
            if len(clmnDescSublist) != 1:
                print("??? INTERNAL ERROR cellTableDescrib - type 'antype' should have no additonal parameter (list)")
                sys.exit(1)
            for atchar in cellTableLine[tcix]:
                if not atchar in antTypeChars:
                    _msgstr = 'Error in input table:\n' + ', '.join(cellTableLine) + '\n'
                    _msgstr += 'Wrong string in column ' + str(tcix+1) + ' (=' + str(cellTableLine[tcix]) + ').'
                    _msgstr += ' Allowed chars: ' + antTypeChars
                    grass.error(_msgstr)
                    err = True

        else:
            print('??? INTERNAL ERROR - invalid type definition in cellTableDescrib ( =', clmnDescSublist[0], ')') 
            sys.exit(1)

if err:
    grass.fatal('Exiting - error(s) found in input table parameters')


# ---- COMPUTATION REGION MANAGEMENT

import subprocess
import atexit


def gen_report_region(regDict):
    """
    Generate region report string:
      n         nn
    w   e  = nn    nn        res = nn    ewres = nn    nres = nn
      s         nn
    """
    nStr = wStr = sStr = eStr = resStr = ewresStr = nsresStr = '?'
    try:
        wStr = str(regDict['w'])
    except:
        pass
    try:
        nStr = str(regDict['n'])
    except:
        pass
    try:
        sStr = str(regDict['s'])
    except:
        pass
    try:
        eStr = str(regDict['e'])
    except:
        pass
    try:
        resStr = str(regDict['res'])
    except:
        pass
    try:
        ewresStr = str(regDict['ewres'])
    except:
        pass
    try:
        nsresStr = str(regDict['nsres'])
    except:
        pass

    lenWstr = len(wStr)
    lenNSstr = max(len(nStr), len(sStr))
    lenEstr = len(eStr)

    msgLine1 = 'N = ' + nStr
    msgLine2 = 'W  E = ' + wStr + ' ' + eStr
    msgLine3 = 'S = ' + sStr
    msgLine4 = 'res E-W/N-S = ' + ewresStr + '/' + nsresStr

    for resKey in ['res', 'ewres', 'nsres']:
        try:
           res = regDict[resKey]
           msgLine2 += resKey + ' = ' + res + '    '
        except:
            pass

    return msgLine1 + '\n' + msgLine2.rstrip() + '\n' + msgLine3 + '\n' + msgLine4


def set_region(regDict):
    regionCmd = ['g.region']
    for key in ['e', 'n', 's', 'w', 'res', 'ewres', 'nsres']:
        try:
            val = regDict[key]
            regionCmd.append(key + '=' + str(val))
        except:
            pass
    if subprocess.call(regionCmd):
        sys.exit(1)


def init_getCSVparam(keyList, paramList):
    """
    Retrieve parameter value (as string) by its header name (string) from the input (CSV) table
    """
    def getCSVparam(key, *args):
        try:
            retstr = paramList[refCellTableHeader.index(key)]
        except (ValueError, IndexError):
            if len(args) == 1:
                if args[0] == 'IGNORE_ERROR':
                    return ''
            print('??? INTERNAL ERROR - unknown key string (parameter name) (=', key, ') while building commands')
            sys.exit(1)
        return retstr
    return getCSVparam


def is_number(s):
     try:
         val=float(s)
         if type(val) is str:  # catch 'nan'
             return False
         return True
     except ValueError:
         return False 

#----------------

# get and save the current region parameters for later restoration
regionOld = grass.region()
atexit.register(set_region, regionOld)

#demFilename = 'dem_ljutomer_25@PERMANENT'
demFilename = options['dem_map']

# create a new region as defined by the region parameters (computation region), and store it
regionStr = options['region']
regionStr = regionStr.replace(':', '=').replace(',', ' ')
if regionStr == 'current':
    pass
else:
    if regionStr == 'dem':
        regionStr = 'raster=' + demFilename
    keyValList = regionStr.split()
    regionCmd = ['g.region']
    for keyVal in keyValList:
        if len(keyVal.split('=')) != 2:
            grass.fatal('Wrong region parameter definition: ' + keyVal)
        [key, val] = keyVal.split('=')
        if key not in ['region', 'raster', 'n', 'e', 's', 'w', 'res']:
            grass.fatal('Unknown region parameter: ' + keyVal)
        if key in ['n', 'e', 's', 'w', 'res']:
            # region is defined numerically
            if not is_number(val):
                grass.fatal('Not a valid parameter value (should be a number): ' + keyVal)
        regionCmd.append(key + '=' + val)
    if subprocess.call(regionCmd):
        grass.fatal('Region parameters: ' + str(regionCmd[2:]))


# "refine" the region: align to resolution, resolution rounded to integer value, the same in both directions

import math

def refineRegCoord(coord, res, mode):
    if not mode in ['+', '-']:
        print("??? INTERNAL ERROR - refineRegCoord mode error (should be '+' or '-')")

    if mode == '+':           # for lower border (W and S)
        x = coord + res/2.0     # x = pixel center (was outer pixel edge (= region border))
        x = math.ceil( x/res) * res  # rounded inwards
        x = x - res/2.0         # back to pixel outer edge

    else:                      # for higher border (N and E)
        x = coord - res/2.0      # x = pixel center (was outer pixel edge (= region border))
        x = math.floor( x/res) * res  # rounded inwards
        x = x + res/2.0          # back to pixel outer edge

    return x

        
regionNotRefined = grass.region()
res = int(round(float(regionNotRefined['ewres'])))
wFl = float(regionNotRefined['w'])
eFl = float(regionNotRefined['e'])
nFl = float(regionNotRefined['n'])
sFl = float(regionNotRefined['s'])
wCompute = refineRegCoord(wFl, res, '+')
eCompute = refineRegCoord(eFl, res, '-')
nCompute = refineRegCoord(nFl, res, '-')
sCompute = refineRegCoord(sFl, res, '+')
regionCmd = ['g.region', 'w=' + str(wCompute), 'e=' + str(eCompute),
             'n=' + str(nCompute), 's=' + str(sCompute), 'res=' + str(res)]
if subprocess.call(regionCmd):
    sys.exit(1)
regionCompute = grass.region()
set_region(regionOld)


# ---- FROM cellTable KEEP ONLY THE CELLS HAVING EFFECT INSIDE THE COMPUTATION REGION ----
# ---- and define the extended computation region (encompasing all the remaining cells)

cellTableReduced = []
wExtended = wCompute
eExtended = eCompute
nExtended = nCompute
sExtended = sCompute
for row in cellTable:
    gCSVp = init_getCSVparam(refCellTableHeader, row)
    cellRadius = int (1000.0 * float(gCSVp('radius')))
    cellEast = int(gCSVp('antEast'))
    cellNorth = int(gCSVp('antNorth'))
    if (cellEast  < wCompute - cellRadius or
        cellEast  > eCompute + cellRadius or
        cellNorth > nCompute + cellRadius or
        cellNorth < sCompute - cellRadius):
        pass
    else:
        cellTableReduced.append(row)
        if cellEast  < wExtended: wExtended = refineRegCoord(cellEast, res, '-')
        if cellEast  > eExtended: eExtended = refineRegCoord(cellEast, res, '+')
        if cellNorth > nExtended: nExtended = refineRegCoord(cellNorth, res, '+')
        if cellNorth < sExtended: sExtended = refineRegCoord(cellNorth, res, '-')
cellTable = cellTableReduced

#define the extended computation region
regionCmd = ['g.region', 'w=' + str(wExtended), 'e=' + str(eExtended),
             'n=' + str(nExtended), 's=' + str(sExtended), 'res=' + str(res)]
if subprocess.call(regionCmd):
    sys.exit(1)
regionExtCompute = grass.region()
set_region(regionOld)


# ---- REMOVE CELLS WITH TRANSMITTER LOCATION NOT DEFINED IN DEM

#testRegionList = [regionExtCompute, regionCompute, regionOld]
testRegionList = [regionExtCompute]
for testRegion in testRegionList:
    set_region(testRegion)
    cellTableReduced = []
    for row in cellTable:
        gCSVp = init_getCSVparam(refCellTableHeader, row)
        cellEastNorthStr = gCSVp('antEast') + ',' + gCSVp('antNorth')
        whatAnsw = grass.read_command('r.what', map = demFilename, coordinates = cellEastNorthStr)
        # re.sub (below) returns height from r.what output, e.g. '590000|151985||212' -> '212'
        heightStr = re.sub(r'.*?\|.*?\|.*?\|', '', whatAnsw)
        if len(heightStr) == 0:
            grass.fatal("Problem with DEM - r.what could not establish transmitter's height")
        if heightStr[0] == '*':
            grassMsg = ('Cell ' + gCSVp('cellName') + '-' + gCSVp('antID') +
                        ' removed from simulation ' + 
                        '(transmitter location ' + cellEastNorthStr + ' not defined in DEM)')
            if testRegion == regionExtCompute: grassMsg += ' (regionExtCompute)'
            elif testRegion == regionCompute: grassMsg += ' (regionCompute)'
            elif testRegion == regionOld: grassMsg += ' (current region)'
            grass.warning(grassMsg)
        else:
            cellTableReduced.append(row)
    cellTable = cellTableReduced
set_region(regionOld)



# ---- PREPARE GRASS_GIS MODEL PROCESSING (hata, ...) ----

grass.info('\nGETTING THE LIST OF EXISTING MODEL AND SECTOR FILES\n' +
             'FROM PREVIOUS SIMULATION RUNS IN THE CURRENT MAPSET...')

# get the list of existing files in current (users's) mapset
p = grass.pipe_command("g.list", type="raster", mapset=mapset)
existingFilesList = p.communicate()[0].splitlines()


cmFilename = options['clutter_map']
# bmFilename = options ['buildings_map']  # for r.urban model

fRecalc = not flags['r']

rxAntHeightStr = options['rx_ant_height']
float( rxAntHeightStr)

requiredModSecFilesList=[]

modelFilenameList=[]
modelCmdList=[]


def param_vals_str(sep):
    pvstr = ''
    ix = 1
    while True:
        pkey = 'P' + str(ix)
        if not pkey in refCellTableHeader:
            break
        pval = gCSVp(pkey,'IGNORE_ERROR')
        if  pval == '':
            break
        else:
            pvstr += sep + pval
        ix += 1
    return pvstr


for row in cellTable:
    gCSVp = init_getCSVparam(refCellTableHeader, row)
    modelName = gCSVp('model')
    modelFilename = ('_' + modelName + param_vals_str('_') + '_' +
                     gCSVp('antEast') + '_' + gCSVp('antNorth') + '_' + gCSVp('antHeightAG') + '_' + 
                     gCSVp('radius') + '_' + gCSVp('freq'))
    if not modelFilename in requiredModSecFilesList:
        requiredModSecFilesList.append(modelFilename)
    if  modelFilename not in modelFilenameList:
        # check for file existence
        if fRecalc or not modelFilename in existingFilesList:
            # file does not exist yet (or recalculation of all files required), proceed
            modelFilenameList.append(modelFilename)

            # model hata (GRASS command r.hata)
            if modelName == 'hata':
                modelCmd = ('r.hata input_dem=' + demFilename + ' output=' + modelFilename +
                            ' area_type=' + gCSVp('P1') + 
                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
                            ' rx_ant_height=' + rxAntHeightStr +
                            ' frequency=' + gCSVp('freq') + ' --overwrite')

            # model hataDEM (GRASS command r.hataDEM)
            elif modelName == 'hataDEM':
                if cmFilename == '':
                    grass.fatal('HataDEM model requires a clutter map to be specified')
                modelCmd = ('r.hataDEM input_dem=' + demFilename + ' clutter=' + cmFilename + ' output=' + modelFilename +
                            ' a0=' + gCSVp('P1') + ' a1=' + gCSVp('P2') + ' a2=' + gCSVp('P3') + ' a3=' + gCSVp('P4') + 
                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
                            ' rx_ant_height=' + rxAntHeightStr +
                            ' frequency=' + gCSVp('freq') + ' --overwrite')

            # model cost231 (GRASS command r.cost231)
            elif modelName == 'cost231':
                modelCmd = ('r.cost231 input_dem=' + demFilename + ' output=' + modelFilename +
                            ' area_type=' + gCSVp('P1') + 
                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
                            ' frequency=' + gCSVp('freq') + ' --overwrite')

            # model Walfisch-Ikegami (GRASS command r.waik)
            elif modelName == 'waik':
                modelCmd = ('r.waik input_dem=' + demFilename + ' output=' + modelFilename +
                            ' free_space_loss_correction=' + gCSVp('P1') + ' bs_correction=' + gCSVp('P2') + 
                            ' range_correction=' + gCSVp('P3') + ' street_width_correction=' + gCSVp('P4') + 
                            ' frequency_correction=' + gCSVp('P5') + ' building_height_correction=' + gCSVp('P6') + 
                            ' street_width=' + gCSVp('P7') + ' distance_between_buildings=' + gCSVp('P8') + 
                            ' building_height=' + gCSVp('P9') + ' phi_street=' + gCSVp('P10') + 
                            ' area_type=' + gCSVp('P11') +
                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
                            ' frequency=' + gCSVp('freq') + ' --overwrite')

            # model fspl (GRASS command r.fspl)
            elif modelName == 'fspl':
                modelCmd = ('r.fspl input_dem=' + demFilename + ' output=' + modelFilename +
                            ' loss_exp=' + gCSVp('P1') + ' loss_offset=' + gCSVp('P2') + 
                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
                            ' rx_ant_height=' + rxAntHeightStr +
                            ' frequency=' + gCSVp('freq') + ' --overwrite')

#            # model itm (GRASS command r.itm)
#            elif modelName == 'itm':
#                modelCmd = ('r.itm input_dem=' + demFilename + ' output=' + modelFilename +
#                            ' relperm=' + gCSVp('P1') + ' conductivity=' + gCSVp('P2') + ' surfref=' + gCSVp('P3') + 
#                            ' radclimate=' + gCSVp('P4') + ' polarization=' + gCSVp('P5') + ' situations=' + gCSVp('P6') + 
#                            ' time=' + gCSVp('P7') + 
#                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
#                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
#                            ' frequency=' + gCSVp('freq') + ' --overwrite')

#            # model urban (GRASS command r.urban)
#            elif modelName == 'urban':
#                if cmFilename == '':
#                    grass.fatal('Urban model requires a clutter map to be specified')
#                if bmFilename == '':
#                    grass.fatal('Urban model requires a buildings map to be specified')
#                modelCmd = ('r.urban input_dem=' + demFilename + ' clutter=' + cmFilename + ' buildings=' + bmFilename +
#                            ' output=' + modelFilename +
#                            ' models_to_use=' + gCSVp('P1') + ' screen_separation=' + gCSVp('P2') +
#                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
#                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
#                            ' frequency=' + gCSVp('freq') + ' --overwrite')

#            # model ITUR1546-4 (GRASS command r.ITUR1546-4)
#            elif modelName == 'ITUR1546-4':
#                modelCmd = ('r.ITUR1546-4 input_dem=' + demFilename +
#                            ' output=' + modelFilename +
#                            ' area_type=' + gCSVp('P1') + ' time_percentage=' + gCSVp('P2') +
#                            ' coordinate=' + gCSVp('antEast') + ',' + gCSVp('antNorth') +
#                            ' ant_height=' + gCSVp('antHeightAG') + ' radius=' + gCSVp('radius') +
#                            ' rx_ant_height=' + rxAntHeightStr +
#                            ' frequency=' + gCSVp('freq') + ' --overwrite')

            else:
                print('??? INTERNAL ERROR - unknown model (=', modelName, ')')

            modelCmdList.append([modelCmd, gCSVp('cellName') + '-' + gCSVp('antID')])


# ---- PREPARE GRASS_GIS SECTOR PROCESSING ----

sectorFilenameList=[]
sectorCmdList=[]
maxPowerSecList=[]

for row in cellTable:
    gCSVp = init_getCSVparam(refCellTableHeader, row)
    modelName = gCSVp('model')
    modelFilename = ('_' + modelName + param_vals_str('_') + '_' +
                     gCSVp('antEast') + '_' + gCSVp('antNorth') + '_' + gCSVp('antHeightAG') + '_' + 
                     gCSVp('radius') + '_' + gCSVp('freq'))
    sectorFilename = (gCSVp('cellName').replace('_','') + '-' + gCSVp('antID') +
                      modelFilename + '_' + gCSVp('antDirection') + '_' + gCSVp('antElecTilt') + '_' +
                      gCSVp('antMechTilt') + '_' + gCSVp('antType'))
    requiredModSecFilesList.append(sectorFilename)
    if sectorFilename in sectorFilenameList:
        print('??? INTERNAL ERROR - duplicated sector in input table (should have been reported)')
        sys.exit(1)
    # prepare data for MaxPower
    maxPowerSec = (gCSVp('cellName') + ';' + gCSVp('antID') + ';' + sectorFilename + '@' + mapset + ';' +
                   gCSVp('power') + ';' + modelName + param_vals_str(';'))
    maxPowerSecList.append(maxPowerSec) 
    # check for file existence
    if fRecalc or not sectorFilename in existingFilesList:
        # file does not exist yet (or recalculation of all files required), proceed
        antennas = findAntennas( antMapTable, gCSVp('antType'), float( gCSVp('freq')), float(gCSVp('antElecTilt')))
        if len(antennas) == 0:
            grass.fatal('No suitable antenna found (type=' + gCSVp('antType') +
                        ' freq=' + gCSVp('freq') + ' etilt=' + gCSVp('antElecTilt') + ')')
##        if len(antennas) > 1:
##            grass.info('Multiple suitable antennas found, the first one will be used:\n' + str(antennas))
        antMSIpath = amDict[antennas[0]];
        sectorFilenameList.append(sectorFilename)
        sectorCmd = ('r.sector pathloss_raster=' + modelFilename + '@' + mapset +
                     ' input_dem=' + demFilename + ' output=' + sectorFilename +
                     ' east=' + gCSVp('antEast') + ' north=' + gCSVp('antNorth') +
                     ' radius=' + gCSVp('radius') +
                     ' ant_data_file=' + antMSIpath + ' beam_direction=' + gCSVp('antDirection') +
                     ' mech_tilt=' + gCSVp('antMechTilt') + ' height_agl=' + gCSVp('antHeightAG') +
                     ' rx_ant_height=' + rxAntHeightStr + ' --overwrite')

        sectorCmdList.append([sectorCmd, gCSVp('cellName') + '-' + gCSVp('antID')])


# ---- DELETE EXISTING FILES NOT REQUIRED BY THIS SIMULATION RUN (for all defined channel models) ----

grass.info('\nDELETING UNNEEDED MODEL AND SECTOR FILES FROM\n' +
             'PREVIOUS SIMULATION RUNS IN THE CURRENT MAPSET...')

fPurge = flags['p']

if fPurge:
    # make a list of defined models (e.g. ['hata', 'hataDEM'])
    modelList = []
    for clmn in cellTableDescrib:
        if clmn[0] == 'model':
            for mlist in clmn[2:]:
               modelList.extend(mlist)

    # list of those files with filenames starting with a defined channel model name,
    #   possibly preceeded by a 'label' (cellName-antennaID) (no '_' allowed)
    # files will be deleted if not needed in this simulation run
    ndel = 0
    delFileList = []
    for fname in existingFilesList:
        fnameNolabel = fname[fname.find('_'):]  # skip 'label0' - find first '_' and check from there on
        for model in modelList:
            if fname.startswith('_' + model + '_') or fnameNolabel.startswith('_' + model + '_'):
                if not fname in requiredModSecFilesList:
                    delFileList.append(fname)
                    ndel += 1
                break
    if ndel == 0:
        grass.info('Purge: no files deleted')
    else:
        _msgstr = 'Purge: the following file are not needed by this run and will be deleted:\n'
        for fname in delFileList:
            _msgstr += fname + '\n'
##        grass.info(_msgstr)
        print(_msgstr)  # grass_info() can not handle very large strings, script crashes

        for fname in delFileList:
            grass.run_command("g.remove",rast=fname)



# ---- COMMANDS EXECUTION (MODEL, SECTOR) ----

def grass_parcmds( cmd_a, ix, spMax):
    """
    Parallel execution of grass command execution (speedup for multicore processors)
    """ 
    global __spList  # subprocess list - local use but it should persist between procedure invocations

    # check if __spList exists, otherwise create it (first calling)
    try:
        __spList
    except NameError:
        __spList = []

    # loop slowly until:
    #  - a new subprocess can be started (limited by max number of subprocesses - spMax)
    #  - and the number of running subprocesses falls (again) below spMAX (usually spMax-1 on exit)
    while True:
        # check for finished subprocesses, delete them from the subprocess list (__spList)
        for ixsp, [sp, ixold] in enumerate(__spList):
            iret = sp.poll()
            if iret != None:
                # a subprocess has finished
                grass.info('< (' + str(ixold+1) + './_)')

                 # print stdout and stderr
                (stdoutdata, stderrdata) = sp.communicate()
                # remove "%" progress substrings (like '  29%<BS><BS><BS><BS><BS>')
                # test example: re.sub( '\ [\ 0-9]{3,3}%\b{5,5}', '', '  29%\b\b\b\b\b')
                stdoutdata = re.sub( b'\ [\ 0-9]{3,3}%[\b]{5,5}', b'', stdoutdata)
                stderrdata = re.sub( b'\ [\ 0-9]{3,3}%[\b]{5,5}', b'', stderrdata)
                # remove final '\n'(s) if any
                stdoutdata = re.sub( b'\n*\Z', b'', stdoutdata)
                stderrdata = re.sub( b'\n*\Z', b'', stderrdata)
                if len(stdoutdata) > 0:
                    grass.info('(O)' + stdoutdata)
                if len(stderrdata) > 0:
                    grass.info('(E)' + stderrdata)

                if iret != 0:
                    # there was an error -> no further processing, just return
                    return iret
                # subprocess finished successfully - clean up and proceed
                __spList.pop(ixsp)

        # if no room available for a new process-> sleep a bit
        # if room available for a new process -> start it and clear cmd_a
        # if no process to start (cmd_a == []) -> return
        if len(__spList) < spMax:
            if cmd_a != []:
               sp = subprocess.Popen(cmd_a, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
               __spList.append([sp,ix])
               cmd_a = []
            else: return 0
        else: time.sleep(0.010)


def exec_cmds(cmdList, parMax):
    """
    Processes commands from commnad list, returns False in case of an error, True otherwise
    """ 
    for ix, [cmd, secName] in enumerate(cmdList):
        grass.info('> ' + secName + ' (' + str(ix+1) + './' + str(len(cmdList)) + ')\n' + cmd)
        # evaluate GISBASE enviroment variable
        cmd_eval = cmd.replace('$GISBASE',os.getenv('GISBASE'))
        cmd_a = cmd_eval.split()
        if parMax == 0:
            # non-parallel execution
            iret = subprocess.call(cmd_a)
        elif parMax > 0:
            # parallel execution - up to parMax parallel executions
            iret = grass_parcmds(cmd_a, ix, parMax)
        if iret != 0:
            return False

    # for parallel execution: wait to complete
    if parMax > 0:
        if grass_parcmds([], 0, 1) != 0:
            return False

    return True

#----------------


grass.info('\nSTARTING RADIO COVERAGE COMPUTATION...')

import time

fCheck = flags['c']

if modelCmdList != []:
    grass.info('          \n----- PROCESSING MODELS -----')

    # set the extended computation region
    grass.info('Temporarily changing the existing current region:\n' + gen_report_region(regionOld) +
               '\nto the extended computation region:\n' + gen_report_region(regionExtCompute))
    set_region(regionExtCompute)
    timeModels0 = time.time()
    if not fCheck:
        if not exec_cmds(modelCmdList, parNum):
            grass.fatal("Error(s) while processing 'model' commands, exiting")
    else:
        grassMsg = 'This is a test run only, the following model commands would be executed:\n'
        for row in modelCmdList:
           grassMsg += row[0] + '\n'
########        grass.info(grassMsg)
        print(grassMsg)  ######## grass.info() seems to fail for large messages (nothing printed)
    timeModels1 = time.time()
    dTimeModels = timeModels1 - timeModels0

# change raster map color table to "rainbow" (which was default in GRASS 6, but not in GRASS 7)
    for rasterMapName in modelFilenameList:
        grass.run_command( "r.colors", map = rasterMapName, color = "rainbow")


if sectorCmdList != []:
    grass.info('          \n----- PROCESSING SECTORS -----')
    timeSectors0 = time.time()
    if not fCheck:
        if not exec_cmds(sectorCmdList, parNum):
            grass.fatal("Error(s) while processing 'sector' commands, exiting")
    else:
        grassMsg = 'This is a test run only, the following model commands would be executed:\n'
        for row in sectorCmdList:
           grassMsg += row[0] + '\n'
########        grass.info(grassMsg)
        print(grassMsg)  ######## grass.info() seems to fail for large messages (nothing printed)
    timeSectors1 = time.time()
    dTimeSectors = timeSectors1 - timeSectors0

# change raster map color table to "rainbow" (which was default in GRASS 6, but not in GRASS 7)
    for rasterMapName in sectorFilenameList:
        grass.run_command( "r.colors", map = rasterMapName, color = "rainbow")


# ---- MAX POWER ----

if maxPowerSecList == []:
    grass.info('Nothing to calculate (no sectors specified)')
else:

    import tempfile

    strOverwrite = ''
##     if grass.gisenv()['OVERWRITE'] == '1': ## supposed to work according to the manual, but it does not
##        strOverwrite = ' --overwrite'
    if "GRASS_OVERWRITE" in os.environ:
        if os.environ["GRASS_OVERWRITE"] == '1':
            strOverwrite = ' --overwrite'

    cellNum = options['cellnum']

    grass.info('          \n----- GENERATING FINAL RESULTS - RASTER MAP AND DB (OPTIONALLY) -----')

    # set the computation region
    grass.info('Temporarily changing current region to the computation region:\n' +
               gen_report_region(regionCompute))
    set_region(regionCompute)

    # Create intermediate (temporary) radio cell/sector table file
    tmpFile = tempfile.NamedTemporaryFile()
    for row in maxPowerSecList:
        tmpFile.write( bytes( row + '\n', 'ascii'))
    tmpFile.flush()
    

    dbDriverName = options['db_driver']
    databaseName = options['database']
    if dbDriverName in ['mysql','pg']:
        dbperfStr = options['dbperf']
    else:
        dbperfStr = '1';

    mpcmd = ("r.MaxPower 'cell_input=" + tmpFile.name + "' output=" + outFilename + 
             ' generate=' + options['generate'] + ' bandwidth=' + options['bandwidth'] +
             ' table=' + dbFilename + ' driver=' + dbDriverName + " 'database=" + databaseName + "' dbperf=" + dbperfStr +
             ' cell_num=' + str(cellNum) + strOverwrite)

    rxThresholdStr = options['rx_threshold']
    if rxThresholdStr != '':
      mpcmd += (' rx_threshold=' + str( rxThresholdStr))


    if dbDriverName == 'none':
        grass.info('> WRITE RASTER MAP ONLY\n' + mpcmd)
    else:
        grass.info('> WRITE RASTER MAP AND DATA TABLE\n' + mpcmd)
    timeWriteDB0 = time.time()
    if fCheck:
        grass.info('This is a test run only, the above r.MaxPower command will not be executed')
    else:
        mpcmd_a = mpcmd.split()
        for ix, param in enumerate(mpcmd_a):
            mpcmd_a[ix] = mpcmd_a[ix].strip("'")
        iret = subprocess.call(mpcmd_a)
        if iret != 0:
            grass.fatal('Error while creating final output files (DB and raster), exiting')
        else:
            # change raster map color table to "rainbow" (which was default in GRASS 6, but not in GRASS 7)
            grass.run_command( "r.colors", map = outFilename, color = "rainbow")
    timeWriteDB1 = time.time()
    dTimeWriteDB = timeWriteDB1 - timeWriteDB0


    # Debugging: display the temporary file
    if fCheck:
        grassMsg = 'A temporary file (' + tmpFile.name + ') has been created for r.MaxPower with the following contents:\n'
        tmpFile.seek(0)
        for row in tmpFile:
            grassMsg += row
##        tmpFile.seek(0)
########        grass.info(grassMsg)
        print(grassMsg)  ######## grass.info() seems to fail for large messages (nothing printed)

##    tmpFile.close()


##grass.info('Restoring the original current region:\n' + gen_report_region(regionOld)
##set_region(regionOld)

if modelCmdList != []:
    grass.info('Processing of %i model(s) took %.2f[s], i.e. %.2f[s/model]'
                % (len(modelCmdList), dTimeModels, dTimeModels/len(modelCmdList)))
if sectorCmdList != []:
    grass.info('Processing of %i sector(s) took %.2f[s], i.e. %.2f[s/sector]'
                % (len(sectorCmdList), dTimeSectors, dTimeSectors/len(sectorCmdList)))
if maxPowerSecList != []:
    if dbDriverName != 'none':
        grass.info('Writing of output raster map and data table took %.2f[s]' % dTimeWriteDB)
    else:
        grass.info('Writing of output raster map took %.2f[s]' % dTimeWriteDB)



# ---- XYZ FILES --- 


if flags['x']:
    grass.info('Writing xyz files...')
    mapNameList = modelFilenameList + sectorFilenameList + [outFilename]
    for mapName in mapNameList:
        fileName = mapName + '.xyz'
        grass.info( ' > ' + fileName)
        cmd = 'r.out.xyz input=' + mapName + ' output=' + fileName
        cmd_a = cmd.split()
        iret = subprocess.call(cmd_a)
        if iret != 0:
            grass.fatal('Error while creating xyz files, exiting')



grass.info('Processing finished')

