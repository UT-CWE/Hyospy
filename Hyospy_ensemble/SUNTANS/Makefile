SUNTANSHOME=/home/dongyu/Documents/Hyospy/CODE/SUNTANS/quad-netcdf-debug/main

include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
  MPIOPTS = 
else
  CC = icc
  MPIDEF = -DNOMPI
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif
ifneq ($(NETCDF4HOME),)
  NETCDFINCLUDE=-I$(NETCDF4HOME)/include
  NETCDFDEF = -DUSENETCDF
else
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFSRC=
endif

LD = $(CC) 
CFLAGS =
MATHLIB = -lm

EXEC=sun
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC) $(NETCDFINCLUDE) 
DEFS = $(MPIDEF) $(NETCDFDEF)
NUMPROCS = 20
datadir=rundata
plotdir=plots/validation

all:	data

test:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) ./$(EXEC) -g -vv --datadir=$(datadir)  
	$(MPIHOME)/bin/mpirun $(MPIOPTS) -np $(NUMPROCS) ./$(EXEC) -s -vvv --datadir=$(datadir)  

debug:	data
	$(MPIHOME)/bin/mpirun -np $(NUMPROCS) xterm -hold -e gdb -command=gdbcommands $(EXEC) 
	#$(MPIHOME)/bin/mpirun -np $(NUMPROCS) --debug $(EXEC) -s -vvv --datadir=$(datadir)

#python ~/suntans/SourceForge/python/SUNTANS/joinsun.py -f GalvCoarse_MayJun8009.nc -p $(datadir) -t 168 -n $(NUMPROCS)
	#python scripts/buildMovies.py
	#python scripts/SuntansFieldComparison.py
	#python scripts/PostProcessSuntans.py

join:  sunjoin
	$(SUNTANSHOME)/sunjoin --steps=24 --np=$(NUMPROCS) --datadir=$(datadir)

validation:
	#ncfile=$(shell awk '/outputNetcdfFile/{print $$2}' $(datadir)/suntans.dat)
	python scripts/SuntansFieldComparison.py $(datadir)/GalvCoarse_MayJun2009_00\*nc $(plotdir)
	python scripts/SuntansADCPComparison.py $(datadir)/GalvCoarse_MayJun2009_00\*nc $(plotdir)

sunjoin: 
	make -C $(SUNTANSHOME) sunjoin
	cp $(SUNTANSHOME)/sunjoin ./
	
data:	initialization.o state.o 
	cp initialization.o state.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME)
	mv $(SUNTANSHOME)/$(EXEC) ./

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	make -C $(SUNTANSHOME)

clean:
	rm -f *.o
	rm -f $(SUNTANSHOME)/*.o

clobber: clean
	rm -rf *~ \#*\# PI* $(EXEC) sunjoin gmon.out rundata/*~
