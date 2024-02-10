# fix this relative to include/
# or use absolute path to the GRASS source code
MODULE_TOPDIR = ../../..

PGM = r.MaxPower

LIBES = $(DBMILIB) $(GISLIB) $(RASTERLIB)
DEPENDENCIES = $(DBMIDEP) $(GISDEP) $(RASTERDEP)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd
