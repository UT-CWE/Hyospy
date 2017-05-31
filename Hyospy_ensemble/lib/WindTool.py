# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 14:03:17 2017

@author: dongyu
"""

"""
WindTool downloads and converts wind data from different sources
"""
from datetime import datetime, timedelta
from netCDF4 import MFDataset, Dataset, num2date
import numpy as np
import utm
import os
import string
import math
import shutil
import urllib
import time
import tarfile
from contextlib import closing

from othertime import SecondsSince
from interpXYZ import interpXYZ

import pdb

class NCEP_wind(object):
    """
    NCEP Surface Winds on Regular Grid
    TAMU server: http://barataria.tamu.edu:8080/thredds/catalog/GNAM-hind-reg-24/catalog.html
    """
    
    def __init__(self, starttime, endtime, **kwargs):
        self.__dict__.update(kwargs)
        
        self.starttime = starttime
        self.endtime = endtime
        
        self.getData()
        #self.writeGNOME(outfile)
        
    def getData(self):
        """
        download the data
        """
        
        basedir = 'http://barataria.tamu.edu:8080/thredds/dodsC/GNAM-hind-reg-24/'
        basefile = 'GNAM-hind-reg'
        
        t1 = datetime.strptime(self.starttime, '%Y-%m-%d-%H')
        t2 = datetime.strptime(self.endtime, '%Y-%m-%d-%H')

        dt = t2 - t1
        dates = []
        for i in range(dt.days + 1):
            dates.append(t1 + timedelta(days=i))
            
        filelist = []
        for t in dates:
            filestr = '%s%s-%s-%s-%s-00-24.nc'%(basedir, basefile, datetime.strftime(t,'%y'), \
                                datetime.strftime(t, '%m'), datetime.strftime(t, '%d'))
            filelist.append(filestr)
        
        vn = MFDataset(filelist, 'r')
        #print vn

        lon = vn.variables['lon'][:]
        lat = vn.variables['lat'][:]
        self.t = vn.variables['time']
        timei = num2date(self.t[:], self.t.units)
        print 'Downloading NCEP wind from %s to %s\n'%(datetime.strftime(timei[0],\
            '%Y-%m-%d %H:%M:%S'), datetime.strftime(timei[-1], '%Y-%m-%d %H:%M:%S'))

        self.tt = timei  # datetime variable        
        self.air_u = vn.variables['air_u'][:]
        self.air_v = vn.variables['air_v'][:]

        #self.y = lat.shape[0]
        #self.x = lon.shape[0]
        self.lon, self.lat = np.meshgrid(lon, lat)

#        ##test plotting 
#        import matplotlib.pyplot as plt
#        from mpl_toolkits.basemap import Basemap
#        west = self.lon.min(); east = self.lon.max()
#        south = self.lat.min(); north = self.lat.max()
#        fig = plt.figure(figsize=(10,10))
#        basemap = Basemap(projection='merc',llcrnrlat=south,urcrnrlat=north,\
#                    llcrnrlon=west,urcrnrlon=east, resolution='h')
#        basemap.drawcoastlines()
#        basemap.fillcontinents(color='coral',lake_color='aqua')
#        basemap.drawcountries()
#        basemap.drawstates()  
#        llons, llats=basemap(self.lon,self.lat)    	
#        basemap.plot(llons, llats,'ob',markersize=8.5)	
#        plt.show()
#        pdb.set_trace()

        
    def writeGNOME(self, outfile):
        """
        Save the data to the file that GNOME can read
        """
        
        GNOME_wind(outfile, self.tt, self.lat, self.lon, self.air_u, self.air_v)
        
        
    def writeSUNTANS(self, outfile):
        """
        Save the data to the file that SUNTANS can read
        """
          
        print 'writing NCEP wind to SUNTANS wind file %s !!!\n'%outfile
        SUNTANS_wind(outfile, self.tt, self.lat, self.lon, self.air_u, self.air_v) 
        
        
class TAMU_NCEP_wind(object):
    """
    TAMU-NCEP wind:
        http://seawater.tamu.edu/tglopu/twdb_lc.tar
    """
    
    def __init__(self, subset_wind=True, **kwargs):
        self.__dict__.update(kwargs)
        
        self.subset_wind = subset_wind
        self.download_file()
        self.getData()
        
    def download_file(self):
        """
        download the data
        """   
    
        basedir = os.getcwd()
        wind_dir = basedir+'/DATA/wind'
        
        if os.path.exists(wind_dir):
            shutil.rmtree(wind_dir)
        os.makedirs(wind_dir) 
        
        url = 'http://seawater.tamu.edu/tglopu/twdb_lc.tar'
        path = '%s/twdb_1c.tar'%wind_dir
        try:
            print 'Opening: %s'%url
            data = urllib.urlopen(url).read()
        except:
            raise Exception, 'cannot open url:\n%s'%url
         
        f = file(path,'wb')
        f.write(data)
        f.close()
        
        # decide whether the download process is over and unzip the tar file
        time.sleep(4)
        if os.path.isfile(path):
            os.mkdir('%s/wind_data'%wind_dir)
            with closing(tarfile.open(path,'r')) as t:
                t.extractall('%s/wind_data'%wind_dir)
        
    def getData(self):
        
        data_dir = os.getcwd()+'/DATA/wind/wind_data'
        
        data=[] # data to store wind data
        #choose the file of the target wind station    
        a=open('%s/twdb000.wndq'%data_dir, 'r').readlines()
        for s in a:
            if '*' not in s:
                if 'days' not in s:
                    line=s.split()
                    data.append(line)

        tt = []        
        for i in range(len(data)):
            tt.append('%s-%s-%s-%s'%(data[i][0], data[i][1], data[i][2], data[i][3]))
            
        self.timei = []
        for t in tt:
            self.timei.append(datetime.strptime(t, '%Y-%m-%d-%H'))            

        print 'Downloading NCEP wind from %s to %s\n'%(datetime.strftime(self.timei[0],\
            '%Y-%m-%d %H:%M:%S'), datetime.strftime(self.timei[-1], '%Y-%m-%d %H:%M:%S'))
                    
        windID, self.lat, self.lon = self.windStations()
        self.air_u = np.zeros((len(self.timei), len(windID)))
        self.air_v = np.zeros((len(self.timei), len(windID)))
            
        for i in range(len(windID)):
            data=[]
            filename = "twdb%03d.wndq"%windID[i]
            a = open('%s/%s'%(data_dir,filename), 'r').readlines()
            for s in a:
                if '*' not in s:
                    if 'days' not in s:
                        line=s.split()
                        data.append(line)
                        for nn in range(len(data)):
                            amp=string.atof(data[nn][4]); angle=string.atof(data[nn][5])
                            self.air_u[nn,i] = amp*math.sin(angle*math.pi/180)
                            self.air_v[nn,i] = amp*math.cos(angle*math.pi/180)            
        
        
    def windStations(self):
        """
        choose the wind file in the right domain
        """
        if self.subset_wind:
            bbox = [-95.3, -94.33, 28.39, 29.88]   # SUNTANS bbox
        else:
            bbox = [-96.03, -93.758, 27.86,29.91]  # ROMS bbox
        
        work_dir=os.getcwd()
        DIR = work_dir+'/DATA/wind/wind_data'
        nfile=len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]) -4
        windID=[]
        lat=[]
        lon=[]
	
        for i in range(nfile):
            filename = "twdb%03d.wndq"%i
            a = open('%s/%s'%(DIR,filename), 'r').readlines()
            for s in a:
                if 'station' in s:
                    line = s.split()
                    lat_tem = string.atof(line[-2])
                    lon_tem = string.atof(line[-1])
                    if bbox[2] < lat_tem < bbox[3] and bbox[0] < lon_tem < bbox[1]:
                        windID.append(i)
                        lat.append(lat_tem)
                        lon.append(lon_tem)
        
        """        
        ######## Plotting the wind station on the basemap figure #########
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        fig = plt.figure(figsize=(10,10))
        basemap = Basemap(projection='merc',llcrnrlat=bbox[2],urcrnrlat=bbox[3],\
                   llcrnrlon=bbox[0],urcrnrlon=bbox[1], resolution='h')
             
        basemap.drawcoastlines()
        basemap.fillcontinents(color='coral',lake_color='aqua')
        basemap.drawcountries()
        basemap.drawstates()  
        
        llons, llats=basemap(lon,lat)
        basemap.plot(llons, llats,'or',markersize=4.5)	
        plt.show() 
        pdb.set_trace()
        """
                                                
        self.bbox = bbox
        return windID, np.asarray(lat), np.asarray(lon)

        
    def writeGNOME(self, outfile):
        """
        Save the data to the file that GNOME can read
        """
            
        ## Step One: intepolate wind to a regular grid
        Num = 20
        lon = np.linspace(self.bbox[0], self.bbox[1], Num)
        lat = np.linspace(self.bbox[2], self.bbox[3], Num) 
        
        lon_new, lat_new = np.meshgrid(lon, lat)
        
        x = np.zeros_like(lat_new)
        y = np.zeros_like(lon_new)
        for i in range(Num):
            for j in range(Num):
                (x[i,j], y[i,j]) = utm.from_latlon(lat_new[i,j], lon_new[i,j])[0:2]
            
        xncep = np.zeros_like(self.lat)
        yncep = np.zeros_like(self.lon)
        for i in range(len(self.lat)):
            (xncep[i], yncep[i]) = utm.from_latlon(self.lat[i], self.lon[i])[0:2]
            
        xy = np.vstack((x.ravel(),y.ravel())).T
        xy_ncep = np.vstack((xncep.ravel(),yncep.ravel())).T
        
        F = interpXYZ(xy_ncep, xy, method='idw')
        
        Nt = len(self.timei)
        air_u_new = np.zeros([Nt, Num, Num])
        air_v_new = np.zeros([Nt, Num, Num])
        
        for tstep in range(Nt):
            utem = F(self.air_u[tstep,:].ravel())
            vtem = F(self.air_v[tstep,:].ravel())
            air_u_new[tstep,:,:]=utem.reshape(Num, Num)
            air_v_new[tstep,:,:]=vtem.reshape(Num, Num)
            
        ## Step Two: write the data to GNOME file
        GNOME_wind(outfile, self.timei, lat_new, lon_new, air_u_new, air_v_new)
        
                
        
        
    def writeSUNTANS(self, outfile):
        """
        Save the data to the file that SUNTANS can read
        """
        
        print 'writing NCEP wind to SUNTANS wind file %s !!!\n'%outfile
        SUNTANS_wind(outfile, self.timei, self.lat, self.lon, self.air_u, self.air_v) 
        


class GNOME_wind(object):
    """
    class that save the data in the wind file of GNOME
    """
    
    def __init__(self, outfile, t, lat, lon, air_u, air_v,**kwargs):
        
        self.__dict__.update(kwargs)
        
        self.outfile = outfile
        self.t = t
        self.lat = lat
        self.lon = lon
        self.air_u = air_u
        self.air_v = air_v
        
        self.y = lat.shape[0]
        self.x = lat.shape[1]
        
        self.writeNC()
        
    def writeNC(self):
        """
        Save the data to the file that GNOME can read
        """
        
        def create_nc_var(outfile, name, dimensions, attdict, dtype='f8',zlib=False,complevel=0):
            
            nc = Dataset(outfile, 'a')
            tmp=nc.createVariable(name, dtype, dimensions,zlib=zlib,complevel=complevel)
            for aa in attdict.keys():
                tmp.setncattr(aa,attdict[aa])
            #nc.variables[name][:] = var	
            nc.close()
        
        #### Write netcdf file ####
        ####create netcdf File####
        nc = Dataset(self.outfile, 'w', format='NETCDF3_CLASSIC')
        nc.file_type = 'Full_Grid'
        nc.Conventions = 'COARDS'
        nc.grid_type = 'curvilinear'
        nc.set_fill_off()
        #nc.Created = datetime.now().isoformat()
        ####Create dimensions####
        nc.createDimension('x', self.x)
        nc.createDimension('y', self.y)
        nc.createDimension('time', None)   ##unlimited dimension
        nc.close()
        
        ####adding variables####
        create_nc_var(self.outfile,'time',('time'),{'units':'seconds since 1970-01-01 00:00:00'})
        create_nc_var(self.outfile,'lon',('y','x'),{'long_name':'Longitude',\
            'units':'degrees_east','standard_name':'longitude'})
        create_nc_var(self.outfile,'lat',('y','x'),{'long_name':'Latitude',\
            'units':'degrees_north','standard_name':'latitude'})
        create_nc_var(self.outfile,'air_u',('time','y','x'),{'long_name':'Eastward Air Velocity',\
            'units':'m/s','missing_value':'-99999.0','standard_name':'eastward_wind'})
        create_nc_var(self.outfile,'air_v',('time','y','x'),{'long_name':'Northward Air Velocity',\
            'units':'m/s','missing_value':'-99999.0','standard_name':'northward_wind'})
        
        time_new = SecondsSince(self.t, basetime = datetime(1970,1,1))
        ######Now writting the variables######
        nc = Dataset(self.outfile,'a')
        nc.variables['time'][:] = time_new
        nc.variables['lon'][:] = self.lon
        nc.variables['lat'][:] = self.lat
        nc.variables['air_u'][:] = self.air_u
        nc.variables['air_v'][:] = self.air_v
        print "Generating NCEP wind file %s for GNOME!!!\n"%self.outfile
        nc.close()
        
             
            
        
class SUNTANS_wind(object):
    """
    class that interpolate, create, and write wind file of SUNTANS
    input: wind data downloaded from other file
    output: convert the data to a file that SUNTANS can read 
    """
    def __init__(self, outfile, t, lat, lon, air_u, air_v,**kwargs):
        
        self.__dict__.update(kwargs)
        
        self.outfile = outfile
        self.t = t
        self.lat = lat
        self.lon = lon
        self.air_u = air_u
        self.air_v = air_v

        ##SUNTANS wind coordinates        
        self.xnarr = np.asarray([ 269237.7775808 ,  299218.17262173,  274549.18690116, 304574.20060701,  \
                        279881.30892322,  309950.57058566,  285233.9569996 ,  315347.0916517 ,  \
                        290606.94066396,  320763.56907416,  296000.06561673,  326199.80428307,  \
                        301413.13371145,  331655.5948561 ,  306845.94294146,  337130.73450553])       
                        
        self.ynarr = np.asarray([ 3100614.94135235,  3095313.75516835,  3130470.64194151,  3125148.90016207, \
                        3160371.25184043,  3155029.14132806,  3190316.40726896,  3184954.11867941, \
                        3220305.7401414 ,  3214923.46793745,  3250338.87805791,  3244936.82052266, \
                        3280415.44429629,  3274993.80354592,  3310535.05780399,  3305094.03980004])
        
        ## Step One: interpolate wind data
        self.interp()
        self.writeNC(self.outfile, self.t, self.xnarr, self.ynarr, self.air_u_new, self.air_v_new)        
        
        
        
    def interp(self):
        """
        interp the data to SUNTANS 16 wind stations
        """
        
        xw=np.zeros_like(self.lat)
        yw=np.zeros_like(self.lon)
        
        if self.lon.ndim == 2:        
            for i in range(self.lon.shape[0]):
                for j in range(self.lon.shape[1]):
                    (xw[i,j],yw[i,j])=utm.from_latlon(self.lat[i,j],self.lon[i,j])[0:2]
        else:
            for i in range(self.lon.shape[0]):
                (xw[i],yw[i])=utm.from_latlon(self.lat[i],self.lon[i])[0:2]
        
        xy_ncep = np.vstack((xw.ravel(),yw.ravel())).T
        xy_narr = np.vstack((self.xnarr.ravel(),self.ynarr.ravel())).T
        
        F = interpXYZ(xy_ncep, xy_narr, method='idw')
        
        Nt = len(self.t[:])
        self.air_u_new = np.zeros([Nt, self.xnarr.shape[0]])
        self.air_v_new = np.zeros([Nt, self.xnarr.shape[0]])
        
        for tstep in range(Nt):
            utem = F(self.air_u[tstep,:].ravel())
            vtem = F(self.air_v[tstep,:].ravel())
            self.air_u_new[tstep,:]=utem
            self.air_v_new[tstep,:]=vtem
            
            
                
    def writeNC(self, outfile, tt, x, y, Uwind, Vwind):
        """
        SUNTANS required wind file, this function creates the netcdf file
        """        
        Nstation = x.shape[0]    
        Nt = len(tt)
        
        nc = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        nc.Description = 'SUNTANS History file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()
        ####Create dimensions####
        nc.createDimension('NVwind', Nstation)
        nc.createDimension('NTair', Nstation)
        nc.createDimension('Nrain', Nstation)
        nc.createDimension('NUwind', Nstation)
        nc.createDimension('NPair', Nstation)
        nc.createDimension('NRH', Nstation)
        nc.createDimension('Ncloud', Nstation)
        nc.createDimension('nt', Nt)
        nc.close()
        
        def create_nc_var(outfile, name, dimensions, attdict, dtype='f8',zlib=False,complevel=0,fill_value=None):
            
            nc = Dataset(outfile, 'a')
            tmp=nc.createVariable(name, dtype, dimensions,zlib=zlib,complevel=complevel,fill_value=fill_value)
            for aa in attdict.keys():
                tmp.setncattr(aa,attdict[aa])
            #nc.variables[name][:] = var	
            nc.close()
        
        ####adding variables####
        create_nc_var(outfile,'x_Vwind',('NVwind'),{'long_name':'Longitude at Vwind','units':'degrees_north'})
        create_nc_var(outfile,'y_Vwind',('NVwind'),{'long_name':'Latitude at Vwind','units':'degrees_east'})
        create_nc_var(outfile,'z_Vwind',('NVwind'),{'long_name':'Elevation at Vwind','units':'m'})
    
        create_nc_var(outfile,'x_Tair',('NTair'),{'long_name':'Longitude at Tair','units':'degrees_north'})
        create_nc_var(outfile,'y_Tair',('NTair'),{'long_name':'Latitude at Tair','units':'degrees_east'})
        create_nc_var(outfile,'z_Tair',('NTair'),{'long_name':'Elevation at Tair','units':'m'})
    
        create_nc_var(outfile,'x_rain',('Nrain'),{'long_name':'Longitude at rain','units':'degrees_north'})
        create_nc_var(outfile,'y_rain',('Nrain'),{'long_name':'Latitude at rain','units':'degrees_east'})
        create_nc_var(outfile,'z_rain',('Nrain'),{'long_name':'Elevation at rain','units':'m'})
        
        create_nc_var(outfile,'x_Uwind',('NUwind'),{'long_name':'Longitude at Uwind','units':'degrees_north'})
        create_nc_var(outfile,'y_Uwind',('NUwind'),{'long_name':'Latitude at Uwind','units':'degrees_east'})
        create_nc_var(outfile,'z_Uwind',('NUwind'),{'long_name':'Elevation at Uwind','units':'m'})
    
        create_nc_var(outfile,'x_Pair',('NPair'),{'long_name':'Longitude at Pair','units':'degrees_north'})
        create_nc_var(outfile,'y_Pair',('NPair'),{'long_name':'Latitude at Pair','units':'degrees_east'})
        create_nc_var(outfile,'z_Pair',('NPair'),{'long_name':'Elevation at Pair','units':'m'})
    
        create_nc_var(outfile,'x_RH',('NRH'),{'long_name':'Longitude at RH','units':'degrees_north'})
        create_nc_var(outfile,'y_RH',('NRH'),{'long_name':'Latitude at RH','units':'degrees_east'})
        create_nc_var(outfile,'z_RH',('NRH'),{'long_name':'Elevation at RH','units':'m'})
    
        create_nc_var(outfile,'x_cloud',('Ncloud'),{'long_name':'Longitude at cloud','units':'degrees_north'})
        create_nc_var(outfile,'y_cloud',('Ncloud'),{'long_name':'Latitude at cloud','units':'degrees_east'})
        create_nc_var(outfile,'z_cloud',('Ncloud'),{'long_name':'Elevation at cloud','units':'m'})
    
        create_nc_var(outfile,'Time',('nt'),{'units':'seconds since 1990-01-01 00:00:00','long_name':'time'})
        create_nc_var(outfile,'Vwind',('nt','NVwind'),{'units':'m s-1','long_name':'Northward wind velocity component','coordinates':'x_Vwind,y_Vwind'})
        create_nc_var(outfile,'Tair',('nt','NTair'),{'units':'Celsius','long_name':'Air Temperature','coordinates':'x_Tair,y_Tair'})
        create_nc_var(outfile,'rain',('nt','Nrain'),{'units':'kg m2 s-1','long_name':'rain fall rate','coordinates':'x_rain,y_rain'})
        create_nc_var(outfile,'Uwind',('nt','NUwind'),{'long_name':'Eastward wind velocity component','coordinates':'x_Uwind,y_Uwind','units':'m s-1'})
        create_nc_var(outfile,'Pair',('nt','NPair'),{'units':'hPa','long_name':'Air Pressure','coordinates':'x_Pair,y_Pair'})
        create_nc_var(outfile,'RH',('nt','NRH'),{'units':'percent','long_name':'Relative Humidity','coordinates':'x_RH,y_RH'})
        create_nc_var(outfile,'cloud',('nt','Ncloud'),{'units':'dimensionless','long_name':'Cloud cover fraction','coordinates':'x_cloud,y_cloud'})
        
            
        z = np.ones([Nstation])*2
        ## change time units
        time_new = SecondsSince(tt)
        ##Tair, rain, Pair, RH, cloud are set to be constant due to a lack of information
        Tair = np.ones([Nt, Nstation])*30.0
        rain = np.ones([Nt, Nstation])*0.0
        Pair = np.ones([Nt, Nstation])*1010.0
        RH = np.ones([Nt, Nstation])*50.0
        cloud = np.ones([Nt, Nstation])*0.0
        ######Now writting the variables######
        nc = Dataset(outfile,'a')
        nc.variables['x_Vwind'][:] = x
        nc.variables['y_Vwind'][:] = y
        nc.variables['z_Vwind'][:] = z
    	
        nc.variables['x_Tair'][:] = x
        nc.variables['y_Tair'][:] = y
        nc.variables['z_Tair'][:] = z
    
        nc.variables['x_rain'][:] = x
        nc.variables['y_rain'][:] = y
        nc.variables['z_rain'][:] = z	
    	
        nc.variables['x_Uwind'][:] = x
        nc.variables['y_Uwind'][:] = y
        nc.variables['z_Uwind'][:] = z
    
        nc.variables['x_Pair'][:] = x
        nc.variables['y_Pair'][:] = y
        nc.variables['z_Pair'][:] = z
    
        nc.variables['x_RH'][:] = x
        nc.variables['y_RH'][:] = y
        nc.variables['z_RH'][:] = z
    
        nc.variables['x_cloud'][:] = x
        nc.variables['y_cloud'][:] = y
        nc.variables['z_cloud'][:] = z
    
        nc.variables['Time'][:] = time_new
        nc.variables['Vwind'][:] = Vwind
        nc.variables['Tair'][:] = Tair
        nc.variables['rain'][:] = rain
        nc.variables['Uwind'][:] = Uwind
        nc.variables['Pair'][:] = Pair
        nc.variables['RH'][:] = RH
        nc.variables['cloud'][:] = cloud
    
        print "Ending writing variables into netcdf file !!!"
        nc.close()
    

        
        
        
        