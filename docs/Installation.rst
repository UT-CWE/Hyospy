Installation Guide For HyosPy
=====

# Documents the packages and libraries required for Hyospy

# UT-Austin CWE

Linux (Tested in 64-bit, Redhat, Ubuntu):

The linux 64bit-python2.7 is required. 

For linux use appropriate package manager (apt-get on Ubuntu, yum on Redhat) to install libraries and dependencies. Pip install is recommended for installing Python packages.


Install SUNTANS:
=====

The first point is to download SUNTANS source code from SUNTANS repository at sourceforge: https://sourceforge.net/projects/suntans/?source=directory. 
The version quad-netcdf-debug is required. More details are in their user guide. 


**Binary Dependencies:**

| gcc, g++ complier
| cmake
| MPI (http://rpm.pbone.net/index.php3/stat/3/srodzaj/1/search/mpich(x86-64)) 
| ParMetis 3.0 (http://glaros.dtc.umn.edu/gkhome/fsroot/sw/parmetis/OLD)
| Triangle (https://www.cs.cmu.edu/~quake/triangle.html)(libX11-devel required)
| NetCDF 4.3.2 or higher (https://github.com/Unidata/netcdf-c/releases/tag/v4.3.2 )(curl, zlib, hdf5 are required to install netcdf library)
| gdal 1.11.1 or higher
| Snappy C library (snappy, snappy-revel required)

Note: gdal library will be installed into "/usr/local/lib" by default. If your system is not yet set up to find libraries in "/usr/local/lib", then you need to add this line to your bashrc (~/.bashrc):

``export LD_LIBRARY_PATH=/usr/local/lib``

Or you can set the environment variable in your current terminal by typing the command in the terminal. Other libraries can be set in a similar way.   


**Python Libraries:**

| setuptools>=23.0
| Numpy=1.11.0
| Scipy
| Cython=0.24.1
| netCDF4-python>=1.0.7
| Ipython
| Matplotlib (libpng-devel, tinter required)
| virtualenv
| basemap>=1.0.7 (geos-3.3.3 required)
| wxPython (streamer required)
| numexpr 


Install pyGNOME:
=====

The linux version of GNOME (General NOAA Operational Modeling Environment) is the oil spill model in Hyospy. The newest version can be found in NOAA repository on Github:
https://github.com/NOAA-ORR-ERD/PyGnome

To make pyGNOME compatible with SUNTANS python libraries, the normal install is preferred. More details are given in their manual, NormalInstall.rst. 


**Python Libraries:**

| pytest
| geojson>=1.3
| gsw=3.0.3
| repoze.lru>=0.6
| pyugrid=0.2.3
| pyshp=1.2
| psutil>=4.3
| pyzmq>=16
| progressbar>=2.3
| testfixtures



*NOAA maintained packages:*

| unit_conversion=2.5.5 (https://github.com/NOAA-ORR-ERD/PyNUCOS)
| py_gd>=0.1.5 (https://github.com/NOAA-ORR-ERD/py_gd.git )(libgd>=2.1.1, https://bitbucket.org/libgd/gd-libgd/downloads)
| pysgrid=0.3.5 (https://github.com/NOAA-ORR-ERD/pysgrid )(cell_tree2d=0.3.0 required, https://github.com/NOAA-ORR-ERD/cell_tree2d.git)


**Oil Library:**

Oil library is required by GNOME, the newest version is available on Github. (https://github.com/NOAA-ORR-ERD/OilLibrary.git)

*Python Libraries:*

| sqlalchemy>=0.7.6
| zope.sqlalchemy>=0.7.6
| awesome-slugify>=1.6
| unidecode>=0.04.19

Note: these packages are required to run Hyospy simulation. More advanced feature of GNOME may require more libraries. Please refer to NOAA website for more details. 


Install TracPy:
=====

TracPy is a second Lagrangian transport model in Hyospy. The source code is available on Github. (https://github.com/kthyng/tracpy)

**Install octant:**

OCTANT (Ocean Modeling Setup and Analysis Tools) is required by TracPy. The newest version is on Github. (https://github.com/hetland/octant)

*Python Libraries:*

pyproj


Install Other dependencies by Hyospy
=====

*Python Libraries*

| utm
| shapely
| pydap
| simplekml

*Binary Dependencies*

netcdf-cxx (ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx-4.2.tar.gz)


**Test**

There are a suite of tests available in GNOME and TracPy. It is recommended to have all the tests passed for the use of Hyospy. 



Install Hyospy
=====

After the requirements are installed, download the Hyospy code from Github. ()

The Python libraries for Hyospy are initiated either in the current directory or subdirectories. So there needs no futher install. Just change a few directories and script mode. 
In subdirectory Hyospy/Hyospy_ensemble/SUNTANS/, change the bash script buildFolder to executable. 
| ``chmod 755 buildFolder``

In file Hyospy/Hyospy_ensemble/SUNTANS/Makefile, modify the directory for SUNTANSHOME. The directory should be where the SUNTANS source code is.

In file Hyospy/Hyospy_ensemble/DATA/data_C/Makefile.in, modify the NETCDF4HOME. The directory is where netcdf library is installed. 
 












