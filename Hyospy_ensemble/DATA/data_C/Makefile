include Makefile.in

OPTFLAGS = -O2

ifneq ($(NETCDF4HOME),)
  NETCDFLIBDIR=-L$(NETCDF4HOME)/lib
  NETCDFINCLUDE=-I$(NETCDF4HOME)/include
  NETCDFLD=-lnetcdf
  NETCDFDEF = -DUSENETCDF
  NETCDFSRC= mynetcdf.c
else
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFDEF = -DNONETCDF
  NETCDFSRC= mynetcdf-nonetcdf.c
endif

#LIBS = $(NETCDFLD)
LIBDIR = $(NETCDFLIBDIR)
#LDFLAGS = -lm $(LIBDIR) $(LIBS)
INCLUDES = $(NETCDFINCLUDE) 
DEFINES = $(NETCDFDEF)
CFLAGS = $(INCLUDES) #$(DEFINES)

LIBS = -lnetcdf_c++ -lnetcdf -lm #-o pgm

CC = g++

all: readnetcdf.cpp
	${CC} -o readnetcdf readnetcdf.cpp $(CFLAGS) $(LIBDIR) $(LIBS)

clean:
	rm -f *.o
