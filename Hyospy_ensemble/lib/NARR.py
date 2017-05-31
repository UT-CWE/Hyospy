# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:42:29 2015

Tools for downloading NARR wind data, create the netcdf file for SUNTANS wind forcing
create the wind file for GNOME
16 wind stations

"""

#from getNARR import getNARR
import numpy as np
import utm
import string
import os
from datetime import datetime, timedelta
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from osgeo import osr
from pydap.client import open_url
from interpXYZ import interpXYZ
import pdb

class windNARR(object):
    """
    general class for downloading NARR wind data and create netcdf file
    """
    def __init__(self,tstart,tend,bbox=[-95.40,-94.30,28.0,29.98],**kwargs):
        
        self.__dict__.update(kwargs)
        self.tstart=tstart
        self.tend=tend
        self.bbox=bbox
        
        self.downloadWind()

	#filename=os.getcwd()+'/CoarseTri/rundata/Galveston_NARR.nc'
	filename='Galveston_NARR.nc'
        self.writeNC(filename)

	#### optional function: plot the wind station ####
	#self.windstation()

	#### function: dump wind data into .WND file (GNOME) ####
	wndfile='wind.WND'
	stationid='11'    ##typically from 0~15
	#self.writeWND(wndfile,stationid)

	#### function: dump wind data into NetCDF grid wind file ####
        
	gridWindFile=os.getcwd()+'/DATA/gridWind.nc'
	#self.gridWind(gridWindFile)

        

    def downloadWind(self):
        """
        download wind from NARR website using Matt's module getNARR
        """
        narr=getNARR(self.tstart,self.tend,self.bbox)
        #############Variables: RH, cloud, Vwind, Uwind, Tair, rain, Pair################
        varlist=['Relative_humidity','Total_cloud_cover','v_wind','u_wind','Temperature','Precipitation_rate','Pressure']
        vardict=narr(varlist)

        var=vardict[0]      #####vardict[0] is the variable dictionary
        (lat, lon)=vardict[1]    #####vardict[1] is the lat, lon
        self.timei=vardict[2]
	
        ######convert the data format for time######
      	self.time=self.SecondsSince(self.timei)
 	#pdb.set_trace()
        ######convert 2D coordination to 1D######
        def _convertM(var):
            """
            Input a matrix (either 2D or 3D)
            squeeze this matrix
            """
            ##2D matrix##
            dim = len(var.shape)
            
            if dim==2:
                return var.ravel()
            elif dim==3:
                return var.reshape(var.shape[0],16)        
        
        self.lat=_convertM(lat)     
        self.lon=_convertM(lon) 
        self.z=np.zeros_like(self.lat)+2.0 
        ######convert lonlat to UTM######
        #(utm_x,utm_y)=utm.from_latlon(28.353786, -95.315109)[0:2]
        for i in range(len(self.lat)):
            (self.lat[i],self.lon[i])=utm.from_latlon(self.lat[i],self.lon[i])[0:2]

        #pdb.set_trace()
        RH=var['Relative_humidity']
        cloud=var['Total_cloud_cover']/100.
        Vwind=var['v_wind']
        Uwind=var['u_wind']
        Tair=var['Temperature']-273.15
        rain=var['Precipitation_rate']
        Pair=var['Pressure']/100.
        ######convert 3D variable to 2D######
        self.RH=_convertM(RH)
        self.cloud=_convertM(cloud)
        self.Vwind=_convertM(Vwind)
        self.Uwind=_convertM(Uwind)
        self.Tair=_convertM(Tair)
        self.rain=_convertM(rain)
        self.Pair=_convertM(Pair)

        ######data length######
        self.Nstation=self.RH.shape[1]
	self.Nt=len(self.time)

	#### variables needed for gridWind function ####
	(self.y,self.x)=lat.shape
	#pdb.set_trace()
        
    def writeNC(self, outfile):
        """
        This function is used to create the netcdf file
        """
        print 'under developing'
        
        ####create netcdf File####
        nc = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        nc.Description = 'SUNTANS History file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()
        ####Create dimensions####
        nc.createDimension('NVwind', self.Nstation)
        nc.createDimension('NTair', self.Nstation)
        nc.createDimension('Nrain', self.Nstation)
        nc.createDimension('NUwind', self.Nstation)
        nc.createDimension('NPair', self.Nstation)
	nc.createDimension('NRH', self.Nstation)
	nc.createDimension('Ncloud', self.Nstation)
	nc.createDimension('nt', self.Nt)
	nc.close()

	####adding variables####
        self.create_nc_var(outfile,'x_Vwind',('NVwind'),{'long_name':'Longitude at Vwind','units':'degrees_north'})
        self.create_nc_var(outfile,'y_Vwind',('NVwind'),{'long_name':'Latitude at Vwind','units':'degrees_east'})
        self.create_nc_var(outfile,'z_Vwind',('NVwind'),{'long_name':'Elevation at Vwind','units':'m'})

	self.create_nc_var(outfile,'x_Tair',('NTair'),{'long_name':'Longitude at Tair','units':'degrees_north'})
        self.create_nc_var(outfile,'y_Tair',('NTair'),{'long_name':'Latitude at Tair','units':'degrees_east'})
        self.create_nc_var(outfile,'z_Tair',('NTair'),{'long_name':'Elevation at Tair','units':'m'})

	self.create_nc_var(outfile,'x_rain',('Nrain'),{'long_name':'Longitude at rain','units':'degrees_north'})
        self.create_nc_var(outfile,'y_rain',('Nrain'),{'long_name':'Latitude at rain','units':'degrees_east'})
        self.create_nc_var(outfile,'z_rain',('Nrain'),{'long_name':'Elevation at rain','units':'m'})

	self.create_nc_var(outfile,'x_Uwind',('NUwind'),{'long_name':'Longitude at Uwind','units':'degrees_north'})
        self.create_nc_var(outfile,'y_Uwind',('NUwind'),{'long_name':'Latitude at Uwind','units':'degrees_east'})
        self.create_nc_var(outfile,'z_Uwind',('NUwind'),{'long_name':'Elevation at Uwind','units':'m'})

	self.create_nc_var(outfile,'x_Pair',('NPair'),{'long_name':'Longitude at Pair','units':'degrees_north'})
        self.create_nc_var(outfile,'y_Pair',('NPair'),{'long_name':'Latitude at Pair','units':'degrees_east'})
        self.create_nc_var(outfile,'z_Pair',('NPair'),{'long_name':'Elevation at Pair','units':'m'})

	self.create_nc_var(outfile,'x_RH',('NRH'),{'long_name':'Longitude at RH','units':'degrees_north'})
        self.create_nc_var(outfile,'y_RH',('NRH'),{'long_name':'Latitude at RH','units':'degrees_east'})
        self.create_nc_var(outfile,'z_RH',('NRH'),{'long_name':'Elevation at RH','units':'m'})

	self.create_nc_var(outfile,'x_cloud',('Ncloud'),{'long_name':'Longitude at cloud','units':'degrees_north'})
        self.create_nc_var(outfile,'y_cloud',('Ncloud'),{'long_name':'Latitude at cloud','units':'degrees_east'})
        self.create_nc_var(outfile,'z_cloud',('Ncloud'),{'long_name':'Elevation at cloud','units':'m'})

	self.create_nc_var(outfile,'Time',('nt'),{'units':'seconds since 1990-01-01 00:00:00','long_name':'time'})
	self.create_nc_var(outfile,'Vwind',('nt','NVwind'),{'units':'m s-1','long_name':'Northward wind velocity component','coordinates':'x_Vwind,y_Vwind'})
	self.create_nc_var(outfile,'Tair',('nt','NTair'),{'units':'Celsius','long_name':'Air Temperature','coordinates':'x_Tair,y_Tair'})
	self.create_nc_var(outfile,'rain',('nt','Nrain'),{'units':'kg m2 s-1','long_name':'rain fall rate','coordinates':'x_rain,y_rain'})
	self.create_nc_var(outfile,'Uwind',('nt','NUwind'),{'long_name':'Eastward wind velocity component','coordinates':'x_Uwind,y_Uwind','units':'m s-1'})
	self.create_nc_var(outfile,'Pair',('nt','NPair'),{'units':'hPa','long_name':'Air Pressure','coordinates':'x_Pair,y_Pair'})
	self.create_nc_var(outfile,'RH',('nt','NRH'),{'units':'percent','long_name':'Relative Humidity','coordinates':'x_RH,y_RH'})
	self.create_nc_var(outfile,'cloud',('nt','Ncloud'),{'units':'dimensionless','long_name':'Cloud cover fraction','coordinates':'x_cloud,y_cloud'})
	

	######Now writting the variables######
	nc = Dataset(outfile,'a')
	nc.variables['x_Vwind'][:]=self.lat
	nc.variables['y_Vwind'][:]=self.lon
	nc.variables['z_Vwind'][:]=self.z
	
	nc.variables['x_Tair'][:]=self.lat
	nc.variables['y_Tair'][:]=self.lon
	nc.variables['z_Tair'][:]=self.z

	nc.variables['x_rain'][:]=self.lat
	nc.variables['y_rain'][:]=self.lon
	nc.variables['z_rain'][:]=self.z	
	
	nc.variables['x_Uwind'][:]=self.lat
	nc.variables['y_Uwind'][:]=self.lon
	nc.variables['z_Uwind'][:]=self.z

	nc.variables['x_Pair'][:]=self.lat
	nc.variables['y_Pair'][:]=self.lon
	nc.variables['z_Pair'][:]=self.z

	nc.variables['x_RH'][:]=self.lat
	nc.variables['y_RH'][:]=self.lon
	nc.variables['z_RH'][:]=self.z

	nc.variables['x_cloud'][:]=self.lat
	nc.variables['y_cloud'][:]=self.lon
	nc.variables['z_cloud'][:]=self.z

	nc.variables['Time'][:]=self.time
	nc.variables['Vwind'][:]=self.Vwind*0.447
	nc.variables['Tair'][:]=self.Tair
	nc.variables['rain'][:]=self.rain
	nc.variables['Uwind'][:]=self.Uwind*0.447
	nc.variables['Pair'][:]=self.Pair
	nc.variables['RH'][:]=self.RH
	nc.variables['cloud'][:]=self.cloud

	print "Ending writing variables into netcdf file !!!"
	nc.close()

	
    

    def create_nc_var(self,outfile, name, dimensions, attdict, dtype='f8',zlib=False,complevel=0,fill_value=None):
        
        nc = Dataset(outfile, 'a')
        tmp=nc.createVariable(name, dtype, dimensions,zlib=zlib,complevel=complevel,fill_value=fill_value)
        for aa in attdict.keys():
            tmp.setncattr(aa,attdict[aa])
        #nc.variables[name][:] = var	
        nc.close()


    def SecondsSince(self,timein,basetime = datetime(1990,1,1)):
        """
        Converts a list or array of datetime object into an array of seconds since "basetime"
    
        Useful for interpolation and storing in netcdf format
        """
        timeout=[]
        try:
            timein = timein.tolist()
        except:
            timein = timein
    
        try:
            for t in timein:
                dt = t - basetime
                timeout.append(dt.total_seconds())
        except:
            dt = timein - basetime
            timeout.append(dt.total_seconds()) 
        
        return np.asarray(timeout)

    def windstation(self):
	"""
	This function is used to plot the wind station from the downloaded NARR data
	"""
	west=-95.40; east=-94.30
	south=27.8;  north=30.0
	fig = plt.figure(figsize=(10,10))
	basemap = Basemap(projection='merc',llcrnrlat=south,urcrnrlat=north,\
		llcrnrlon=west,urcrnrlon=east, resolution='h')
          
	basemap.drawcoastlines()
    	basemap.fillcontinents(color='coral',lake_color='aqua')
    	basemap.drawcountries()
    	basemap.drawstates()  
        
	lat=np.zeros_like(self.lat)
	lon=np.zeros_like(self.lon)
	for i in range(len(self.lat)):
		(lat[i],lon[i])=utm.to_latlon(self.lat[i],self.lon[i],15,'U')
    	llons, llats=basemap(lon,lat)
    	basemap.plot(llons, llats,'or',markersize=4.5)	
	plt.show()


    def writeWND(self,wndfile,stationid):
	"""
	This function is used to dump the wind data into .WND file (needed by GNOME)
	Note: The .WND file is one-value Data set for wind forcing

	wndfile: the output filename
	stationid: the station where the value is used for wind data
	"""
	ID=string.atoi(stationid)
	Uwind=self.Uwind[:,ID]
	Vwind=self.Vwind[:,ID]
	amp=np.sqrt(np.power(Uwind,2)+np.power(Vwind,2))*0.447 ## wind amplitude
	angle=np.degrees(np.arctan2(Vwind,Uwind))  ## wind angle in degrees

	lat=self.lat[ID]
	lon=self.lon[ID]
	(lat,lon)=utm.to_latlon(lat,lon,15,'U')	
	
	nt=len(self.timei)
	data0=np.zeros((nt,7)).tolist()

	for tt in range(nt):
	    ## wind direction
	    data0[tt][6]=str(angle[tt])
	    ## wind amplitude
	    data0[tt][5]=str(amp[tt])
   	    ## time
	    data0[tt][4]='00,'
	    data0[tt][3]=str('%02d' % self.timei[tt].hour)+','
	    data0[tt][2]=str('%02d' % self.timei[tt].year)[2:]+','
	    data0[tt][1]=str('%01d' % self.timei[tt].month)+','
	    data0[tt][0]=str('%01d' % self.timei[tt].day)+','
	

	dire=os.getcwd()
	if os.path.isfile(dire+'/GNOME/wind.WND'):
            print 'The GNOME wind data file already exists, generate new file!!!'
            os.remove(dire+'/GNOME/wind.WND')
        else:
            print 'Generating GNOME wind file!!'
        
        ff=open(dire+'/wind.WND','w')

        ff.write('NARR'+str(ID)+'\n')
        ff.write(str(lat)+' '+str(lon)+'\n')
        ff.write('mps'+'\n')
        ff.write('LTime'+'\n')
        ff.write('0,0,0,0,0,0,0,0'+'\n')
        for i in (data0):
            k='   '.join([str(j) for j in i])
            ff.write(k+'\n')
        ff.close
		
	
    def gridWind(self,outfile):
	"""
	Instead of single point wind, this function is used to created NetCDF curvilinear grid 
	surface wind model
	"""
	#### Prepare data for GNOME gridWind netcdf file ####
	#### utm to degree ####
	lat = np.zeros_like(self.lat)
	lon = np.zeros_like(self.lon)

	for i in range(len(self.lat)):
            (lat[i],lon[i])=utm.to_latlon(self.lat[i],self.lon[i], 15, 'R')

	lon = lon.reshape(self.y, self.x)
	lat = lat.reshape(self.y, self.x)
	air_u = self.Uwind.reshape(len(self.Uwind),self.y,self.x)
	air_v = self.Vwind.reshape(len(self.Vwind),self.y,self.x)
	time =  self.SecondsSince(self.timei, basetime = datetime(1970,1,1))	

	####create netcdf File####
        nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')   #### be careful about the file type
        nc.file_type = 'Full_Grid'
        nc.Conventions = 'COARDS'
	nc.grid_type = 'curvilinear'
	nc.title = 'Forecast:wind+tide+river'
	nc.set_fill_off()
        #nc.Created = datetime.now().isoformat()
        ####Create dimensions####
        nc.createDimension('x', self.x) 
        nc.createDimension('y', self.y)
        nc.createDimension('time', None)   ##unlimited dimension
	nc.close()

	#pdb.set_trace()
	####adding variables####
        self.create_nc_var(outfile,'time',('time'),{'long_name':'Time',\
		'units':'seconds since 1970-01-01 00:00:00'})
	self.create_nc_var(outfile,'lon',('y','x'),{'long_name':'Longitude',\
		'units':'degrees_east','standard_name':'longitude'})
	self.create_nc_var(outfile,'lat',('y','x'),{'long_name':'Latitude',\
		'units':'degrees_north','standard_name':'latitude'})
	self.create_nc_var(outfile,'air_u',('time','y','x'),{'long_name':'Eastward Air Velocity',\
		'units':'m/s','missing_value':'-99999.0','standard_name':'eastward_wind'})
	self.create_nc_var(outfile,'air_v',('time','y','x'),{'long_name':'Northward Air Velocity',\
		'units':'m/s','missing_value':'-99999.0','standard_name':'northward_wind'})

	######Now writting the variables######
	nc = Dataset(outfile,'a')
	nc.variables['time'][:] = time
	nc.variables['lon'][:] = lon
	nc.variables['lat'][:] = lat
	nc.variables['air_u'][:] = air_u
	nc.variables['air_v'][:] = air_v
	print "Ending writing variables into %s netcdf file !!!"%outfile
	nc.close()

	#pdb.set_trace()	

        
class getNARR(object):
    """
    Tool for downloading data from the North American Regional Reanalysis (NARR) product

    Note that this requires the netCDF4 module to be compiled with opendap support (libcurl)

    Created on Thu Nov 08 17:22:40 2012

    @author: mrayson

    """
    """
    General class for initialising a NARR data object
    
    Example Usage:
        tstart = '20100101'
        tend = '20100103'
        bbox = [-95.40,-94.49,28.8,29.9]
        varname = 'Relative_humidity'
        
        narr = getNARR(tstart,tend,bbox)
        RH = narr(varname)
        
    """
    basedir = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/narr/'
    basefile = 'narr-a_221_'
    verbose = True
    
    def __init__(self,tstart,tend,bbox,**kwargs):
        
        self.__dict__.update(kwargs)
        self.tstart=tstart
        self.tend=tend
        self.bbox=bbox
        
        # Get the filenames and the coordinates
        self.getFileNames()
        
        self.getNARRlatlon()
        
        self.findXYindices()
        
        self.getArraySize()
        
    def __call__(self,varnames):
        """
        Call function to extract variables in list, "varnames"
        """

        return self.getDataPyDap(varnames), (self.lat, self.lon), self.time
        #return self.getData(varnames)
    
    def getData(self,varnames):
        """
        Downloads the variable data using the NetCDF4 module
        
        This module has memory leaks if not compiled with the correct netcdf/hdf/curl libraries
        
        Safer to use pydap.
        """
        
        
        data={}
        self.nc = Dataset(self.grbfiles[0],'r')
        for vv in varnames:
            if self.verbose:
                print 'Retrieving variable: %s dimension info...'%vv        
            self.getDimInfo(vv)
            data.update({vv:np.zeros((self.nt,self.ny,self.nx))})
        self.nc.close()
        
        tt=-1
        for ff in self.grbfiles:
            tt+=1
            if self.verbose:
                print '     File: %s...'%ff
                
            nc = Dataset(ff,'r')
            
            for vv in varnames:
                print vv
                # get the height coordinate for 4D arrays
                if self.ndim == 4:
                     # dimension order [time, z, y, x]
                    data[vv][tt,:,:] = nc.variables[vv][:,0,self.y1:self.y2,self.x1:self.x2]
                elif self.ndim == 3:
                    data[vv][tt,:,:] = nc.variables[vv][:,self.y1:self.y2,self.x1:self.x2]
                
            nc.close()    
        
        return data
        
    def getDataPyDap(self,varnames):
        """
        Downloads the variable data using the pydap library
        """
        
        data={}
        dims={}
        self.nc = Dataset(self.grbfiles[0],'r')
        for vv in varnames:
            if self.verbose:
                print 'Retrieving variable: %s dimension info...'%vv  
		     
            self.getDimInfo(vv)
            data.update({vv:np.zeros((self.nt,self.ny,self.nx))})
            dims.update({vv:self.ndim})
        self.nc.close()
        

        tt=-1
        for ff in self.grbfiles:
            tt+=1
            if self.verbose:
                print '     File: %s...'%ff
                
            nc = open_url(ff)

            for vv in varnames:
                print vv
                self.ndim  = dims[vv]
                # get the height coordinate for 4D arrays
                if self.ndim == 4:
                     # dimension order [time, z, y, x]
                    data[vv][tt,:,:] = nc[vv][:,0,self.y1:self.y2,self.x1:self.x2]
                elif self.ndim == 3:
                    data[vv][tt,:,:] = nc[vv][:,self.y1:self.y2,self.x1:self.x2]
                   
        
        return data
        
    def getDimInfo(self,vv):
        """
        Gets the variable dimension data
        """
       
        nc = self.nc
	
        self.dimdata = nc.variables[vv].dimensions
        self.ndim  = nc.variables[vv].ndim
        
        # get the height coordinate for 4D arrays
        if self.ndim == 4:
             # dimension order [time, z, y, x]
            self.z = nc.variables[self.dimdata[1]][0]
        else:
            self.z=0.0
    
    def getArraySize(self):
        
        self.nt = len(self.time)
        self.nx = self.x2-self.x1
        self.ny = self.y2-self.y1
        
    def getFileNames(self):
        """
        Return the filenames for each time step between the two dates
        """

        # build a list of timesteps
        t1 = datetime.strptime(self.tstart,'%Y%m%d')
        t2 = datetime.strptime(self.tend,'%Y%m%d')
        
        self.time = []
        t0=t1
        while t0 < t2:
            self.time.append(t0)
            t0 += timedelta(hours=3)
            
        # Build the list of filenames
        self.grbfiles=[]
        for tt in self.time:
            filestr='%s%s/%s/%s%s_000.grb'%(self.basedir,datetime.strftime(tt,'%Y%m'),\
            datetime.strftime(tt,'%Y%m%d'),self.basefile,datetime.strftime(tt,'%Y%m%d_%H%M'))
            self.grbfiles.append(filestr)
            
    def findXYindices(self):
        """
        Find the indices of the lower left and upper right corners from the grid
        """
        
        dist = np.sqrt( (self.lon - self.bbox[0])**2 + (self.lat - self.bbox[2])**2)
        j = np.argwhere(dist==dist.min())
        self.y1=j[0,0]
        self.x1=j[0,1]
        
        dist = np.sqrt( (self.lon - self.bbox[1])**2 + (self.lat - self.bbox[3])**2)
        j = np.argwhere(dist==dist.min())
        self.y2=j[0,0]
        self.x2=j[0,1]
        
        # Resize the lat lon arrays
        self.lat=self.lat[self.y1:self.y2,self.x1:self.x2]
        self.lon=self.lon[self.y1:self.y2,self.x1:self.x2]
	

    def getNARRlatlon(self):
        """
        Returns the NARR grid in lat/lon coordinates
        
        Need to convert from their Lambert Conformal projection to WGS84
        
        *** THIS NEEDS TO BE CHECKED THOROUGHLY!!! ***
        """
        
        nc = Dataset(self.grbfiles[0],'r')

        x = nc.variables['x'][:]*1000.0
        y = nc.variables['y'][:]*1000.0
        
        nc.close()
        
        # NARR grid projection information:
        # -----------------------------------    
        #grid_mapping_name: lambert_conformal_conic
        #standard_parallel: 50.0
        #longitude_of_central_meridian: -107.0
        #latitude_of_projection_origin: 50.0
        #earth_shape: spherical
        #earth_radius: 6367470.21484375
        #GRIB_param_Dx: 32463.0
        #GRIB_param_Dy: 32463.0
        #GRIB_param_GDSkey: 55295
        #GRIB_param_La1: 1.0
        #GRIB_param_Latin1: 50.0
        #GRIB_param_Latin2: 50.0
        #GRIB_param_Lo1: -145.5
        #GRIB_param_LoV: -107.0
        #GRIB_param_NpProj: true
        #GRIB_param_Nx: 349
        #GRIB_param_Ny: 277
        #GRIB_param_ProjFlag: 0
        #GRIB_param_ResCompFlag: 8
        #GRIB_param_SpLat: 0.0
        #GRIB_param_SpLon: 0.0
        #GRIB_param_VectorComponentFlag: gridRelative
        #GRIB_param_Winds: Relative
        #GRIB_param_grid_name: Lambert_Conformal
        #GRIB_param_grid_radius_spherical_earth: 6367.47
        #GRIB_param_grid_shape: spherical
        #GRIB_param_grid_shape_code: 0
        #GRIB_param_grid_type: 3
        #GRIB_param_grid_units: m
        #GRIB_param_scanning_mode: 64
        #DODS:
        #  strlen: 0
        
        # Define the input coordinate system
        # See this website:
            # http://trac.osgeo.org/proj/wiki/GenParms
        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=lcc +lat_0=50.0 +lon_0=-107.0 +lon_1=-145.5 +lat_1=50.0 +a=6367470.21484375 +b=6367470.21484375 +units=m')    
        
        # set the output coordinate system
        srsout = osr.SpatialReference()
        srsout.SetWellKnownGeogCS( "WGS84" );
        
        # define the transformation object
        ct = osr.CoordinateTransformation(srs, srsout)
        
        nx = len(x)
        ny = len(y)
        self.lon = np.zeros((ny,nx))
        self.lat = np.zeros((ny,nx)) 
        for ii in range(nx):
            for jj in range(ny):
                X,Y,z =  ct.TransformPoint(x[ii],y[jj])
                self.lon[jj,ii]=X
                self.lat[jj,ii]=Y
        
def interp_wind(current_file):
    """
    This function interpolates the wind dataset (from 16 stations for example) into
    the blended grid of current velocity
    """
    base_dir=os.getcwd()
    #current_file=base_dir+'/DATA/blendeduv.nc'
    nc1=Dataset(current_file,'r')
    #print nc1

    lon=nc1.variables['lon'][:]
    lat=nc1.variables['lat'][:]
    (y,x)=lon.shape

    xroms=np.zeros_like(lon)  ## xroms: longitude, yroms: latitude
    yroms=np.zeros_like(lat)
    for i in range(y):
        for j in range(x):
            (yroms[i,j],xroms[i,j])=utm.from_latlon(lat[i,j],lon[i,j])[0:2]
        

    wind_file=base_dir+'/DATA/gridWind.nc'
    nc2=Dataset(wind_file,'r')
    #print nc2

    lonw=nc2.variables['lon'][:]
    latw=nc2.variables['lat'][:]
    time=nc2.variables['time'][:]
    air_u=nc2.variables['air_u'][:]
    air_v=nc2.variables['air_v'][:]
    (y2,x2)=lonw.shape

    xw=np.zeros_like(lonw)  ## xroms: longitude, yroms: latitude
    yw=np.zeros_like(latw)
    for i in range(y2):
        for j in range(x2):
            (yw[i,j],xw[i,j])=utm.from_latlon(latw[i,j],lonw[i,j])[0:2]
        
    xy_roms=np.vstack((xroms.ravel(),yroms.ravel())).T
    xy_wind=np.vstack((xw.ravel(),yw.ravel())).T

    Fuv=interpXYZ(xy_wind,xy_roms,method='idw')

    uout=np.zeros([len(time),y,x])
    vout=np.zeros([len(time),y,x])

    for tstep in range(len(time)):
        utem=Fuv(air_u[tstep,:,:].ravel()).reshape(y,x)
        vtem=Fuv(air_v[tstep,:,:].ravel()).reshape(y,x)
        uout[tstep,:,:]=utem
        vout[tstep,:,:]=vtem

    #plt.pcolor(uout[3,:,:])
    #plt.show()

#############################################################################
#############################################################################
    def create_nc_var(outfile, name, dimensions, attdict, dtype='f8',zlib=False,complevel=0):
            
        nc = Dataset(outfile, 'a')
        tmp=nc.createVariable(name, dtype, dimensions,zlib=zlib,complevel=complevel)
        for aa in attdict.keys():
            tmp.setncattr(aa,attdict[aa])
        #nc.variables[name][:] = var	
        nc.close()
    #### Write netcdf file ####
    ####create netcdf File####
    outfile=base_dir+'/GNOME/interp_wind.nc'
    nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
    nc.file_type = 'Full_Grid'
    nc.Conventions = 'COARDS'
    nc.grid_type = 'curvilinear'
    nc.set_fill_off()
    #nc.Created = datetime.now().isoformat()
    ####Create dimensions####
    nc.createDimension('x', x) 
    nc.createDimension('y', y)
    nc.createDimension('time', None)   ##unlimited dimension
    nc.close()

    #pdb.set_trace()
    ####adding variables####
    create_nc_var(outfile,'time',('time'),{'units':'seconds since 1970-01-01 00:00:00'})
    create_nc_var(outfile,'lon',('y','x'),{'long_name':'Longitude',\
        'units':'degrees_east','standard_name':'longitude'})
    create_nc_var(outfile,'lat',('y','x'),{'long_name':'Latitude',\
        'units':'degrees_north','standard_name':'latitude'})
    create_nc_var(outfile,'air_u',('time','y','x'),{'long_name':'Eastward Air Velocity',\
        'units':'m/s','missing_value':'-99999.0','standard_name':'eastward_wind'})
    create_nc_var(outfile,'air_v',('time','y','x'),{'long_name':'Northward Air Velocity',\
        'units':'m/s','missing_value':'-99999.0','standard_name':'northward_wind'})

    ######Now writting the variables######
    #uout = np.ones_like(uout) #just for testing if grid wind is added
    #vout = np.ones_like(vout)
    nc = Dataset(outfile,'a')
    nc.variables['time'][:] = time
    nc.variables['lon'][:] = lon
    nc.variables['lat'][:] = lat
    nc.variables['air_u'][:] = uout
    nc.variables['air_v'][:] = vout
    print "Ending writing variables into %s netcdf file !!!"%outfile
    nc.close()    
         

#############
# Testing stuff

#tstart = '20100101'
#tend = '20100102'
#bbox = [-95.40,-94.49,28.8,29.9]
##
#
#narr=getNARR(tstart,tend,bbox)
#
##RH = narr('Relative_humidity')
#cloud = narr('Total_cloud_cover')        
        
###########################For Testing############################
if __name__ == "__main__":
    tstart = '20140301'
    tend = '20140330'
    bbox = [-95.40,-94.30,28.0,29.98]   #This bbox generates 16 wind stations
    windNARR(tstart, tend, bbox)

    print "end"





