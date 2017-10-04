# -*- coding: utf-8 -*-
"""
Create the fully coupled Galveston Input files

Created on Mon Apr 01 18:00:32 2013

@author: mrayson
"""

def generateBI(start,end, OBC_opt='file', IC_opt='constant'):

    from sundriver import sundriver, dumpinputs

    ####
    # Variable inputs for the class
    #starttime = '20070301.000000'
    #endtime = '20071001.000000'
    starttime=start
    endtime=end
    dt = 900.0
    sunpath = 'SUNTANS/rundata'

    plotdir = 'SUNTANS/plots/inputs'
    ####

    # Initiate the driver class
    sun = sundriver()

    # Switches to generate bathymetry, boundary, meteorology and initial condition input files
    sun.makebathy=False
    sun.makebnd=True
    sun.makewinds=False
    sun.makeinitial=True

    sun.modifyedges =True# Option to modify the boundary edges

    ####
    # General options
    ####
    # Grid projection variables
    sun.convert2utm=False
    sun.CS='NAD83'
    sun.utmzone=15
    sun.isnorth=True
    sun.vdatum = 'MSL'

    # Verical grid options
    sun.Nkmax = 20 # number of layers
    sun.r = 1.08 # vertical stretching parameter

    ####
    # Bathymetry interpolation options
    ####
    sun.depthfile = 'DATA/NOAA_25m_UTM_DEM.nc'
    sun.depthmax=0.1
    sun.interpmethod='kriging' # Interpolation method:  'nn', 'idw', 'kriging', 'griddata'
    sun.plottype='mpl' # Type of plot: 'mpl', 'vtk2' or 'vtk3'

    # Interpolation options
    sun.NNear=10

    # Interpolate to nodes then take maximum depth for cell
    sun.interpnodes=False

    # IF interpmethod = 'idw' 
    sun.p = 1.0 #  power for inverse distance weighting

    # IF interpmethod = 'kriging' 
    sun.varmodel = 'spherical'
    sun.nugget = 0.1
    sun.sill = 0.8
    sun.vrange = 2500.00

    # Smoothing options
    sun.smooth=True
    sun.smoothmethod='kriging' # USe kriging or idw for smoothing
    sun.smoothnear=6 # No. of points to use for smoothing

    #####
    # Open boundary options
    ####
    sun.opt_bcseg = 'file' # Segment boundary condition option: 'constant' or 'file'
    sun.opt_bctype2 = 'constant' # Type 2 boundary condition option: 'constant'
    #sun.opt_bctype3 = 'ROMSFILE' # Type 3 boundary condition option: 'constant',
    #,'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'
    sun.opt_bctype3 = OBC_opt

    sun.bcpolygonfile = 'SUNTANS/gis/GalvCoarse_BndPoly_Rivs_TWDB.shp' # Shape file with fields 'marker' and 'edge_id'
    sun.bcfile = 'GalvCoarse_BC.nc' # Input boundary condition file
    
    # IF opt_bcseg = 'consant
    sun.Q0 = 5.0 # m3/s

    # IF opt_bctype2/opt_bctype3 = 'constant'
    sun.T0 = 30.0 # Open boundary and initial condition background temperature
    sun.S0 = 0.0 # Open boundary and initial condition background salinity

    # IF opt_bctype2 = 'file' or 'ROMSFILE'
    sun.waterlevelstationID = 8772447 # USCG Freeport gage
    #n.waterlevelstationID = 8771341 # Galveston Entrance North Jetty
    #sun.waterlevelstationID = 8771450 # Galveston Pier 21

    # IF opt_bytype2 = 'filt'
    sun.TairstationID = '722436-12906' # Houston-Ellington weather station

    sun.filttype='low'
    sun.cutoff = 7200.0

    # Option to zero ROMS uv and eta at boundaries and initial conditions
    sun.useROMSuv=True
    sun.useROMSeta=False

    # Option to use OTS velocities along bounaries
    sun.useOTISuv=True
    ####
    # Initial condition options
    ####
    #sun.opt_ic = 'ROMS' #'constant', 'depth_profile', 'ROMS', 'SUNTANS'
    sun.opt_ic = IC_opt

    sun.suntansicfile = None

    sun.agesourcepoly = None

    sun.icfile = 'GalvCoarse_IC.nc'

    sun.icfilterdx = 0.0 # spatial filter length scale for smoothing ic's
    ####
    # Metfile
    #### 

    # For plotting only
    sun.metfile = 'Galveston_NARR_20122016.nc'
    ####
    # Input file names
    #### 
    sun.romsfile = ['DATA/txla_subset_HIS.nc']
    sun.otisfile = 'DATA/Tides/Model_Mex'
    sun.dbasefile = 'DATA/GalvestonObs.db'

    ######################
    # Now call the class...
    ######################

    sun(sunpath,starttime,endtime,dt)


def generateBI2(start,end):
    from sundriver import sundriver, dumpinputs

###
# Variable inputs for the class
    starttime = start
    endtime = end
    dt = 900.0
    sunpath = 'FineTri/rundata'

    plotdir = 'FineTri/plots/inputs'
###

# Initiate the driver class
    sun = sundriver()

# Switches to generate bathymetry, boundary, meteorology and initial condition input files
    sun.makebathy=False
    sun.makebnd=True
    sun.makewinds=False
    sun.makeinitial=True

    sun.modifyedges =True # Option to modify the boundary edges

###
# General options
###
# Grid projection variables
    sun.convert2utm=False
    sun.CS='NAD83'
    sun.utmzone=15
    sun.isnorth=True
    sun.vdatum = 'MSL'

# Verical grid options
    sun.Nkmax = 20 # number of layers
    sun.r = 1.08 # vertical stretching parameter

###
# Bathymetry interpolation options
###
    sun.depthfile = 'DATA/NOAA_25m_UTM_DEM.nc'
    sun.depthmax=0.1
    sun.interpmethod='kriging' # Interpolation method:  'nn', 'idw', 'kriging', 'griddata'
    sun.plottype='mpl' # Type of plot: 'mpl', 'vtk2' or 'vtk3'

# Interpolation options
    sun.NNear=6

# Interpolate to nodes then take maximum depth for cell
    sun.interpnodes=False

# IF interpmethod = 'idw' 
    sun.p = 1.0 #  power for inverse distance weighting

# IF interpmethod = 'kriging' 
    sun.varmodel = 'spherical'
    sun.nugget = 0.1
    sun.sill = 0.8
    sun.vrange = 2500.00

# Smoothing options
    sun.smooth=False
    sun.smoothmethod='kriging' # USe kriging or idw for smoothing
    sun.smoothnear=4 # No. of points to use for smoothing


# option to adjust channel depths using a shapefile
    sun.adjust_depths=True
    sun.channel_shpfile='gis/channel_depths.shp'


####
# Open boundary options
####
    sun.opt_bcseg = 'file' # Segment boundary condition option: 'constant' or 'file'
    sun.opt_bctype2 = 'constant' # Type 2 boundary condition option: 'constant'
    sun.opt_bctype3 = 'file' # Type 3 boundary condition option: 'constant',
#,'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'

    sun.bcpolygonfile = 'FineTri/gis/GalvQuadTri_BndPoly_Rivs_TWDB_tri.shp' # Shape file with fields 'marker' and 'edge_id'
    sun.bcfile = 'GalvTri_BC.nc' # Input boundary condition file
    
# IF opt_bcseg = 'consant
    sun.Q0 = 0.0 # m3/s

# IF opt_bctype2/opt_bctype3 = 'constant'
    sun.T0 = 30.0 # Open boundary and initial condition background temperature
    sun.S0 = 0.0 # Open boundary and initial condition background salinity

# IF opt_bctype2 = 'file' or 'ROMSFILE'
#sun.waterlevelstationID = 8772447 # USCG Freeport gage
    sun.waterlevelstationID = 8771450 # Galveston Pier 21

# IF opt_bytype2 = 'filt'
    sun.TairstationID = '722436-12906' # Houston-Ellington weather station

    sun.filttype='low'
    sun.cutoff = 7200.0

# Option to zero ROMS uv and eta at boundaries and initial conditions
    sun.useROMSuv=False
    sun.useROMSeta=False

# Option to use OTS velocities along bounaries
    sun.useOTISuv=False
####
# Initial condition options
####
    sun.opt_ic = 'ROMS' #'constant', 'depth_profile', 'ROMS'

    sun.icfile = 'GalvTri_IC.nc'

    sun.icfilterdx = None # spatial filter length scale for smoothing ic's

    sun.agesourcepoly = None
#sun.agesourcepoly = 'gis/AgeSourcePoly_Entrance.shp'
###
# Metfile
### 

# For plotting only
    sun.metfile = 'Galveston_NARR_2010.nc'
###
# Input file names
### 
    sun.romsfile = ['DATA/txla_subset_HIS_tide123.nc']
    sun.otisfile = 'DATA/Tides/Model_Mex'
    sun.dbasefile = 'DATA/GalvestonObs.db'

######################
# Now call the class...
######################

    sun(sunpath,starttime,endtime,dt)


####
## Dump figures of the input data
#dump = dumpinputs(suntanspath=sunpath,icfile=sun.icfile,bcfile=sun.bcfile,metfile=sun.metfile)

#dump(plotdir)

