# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:24:26 2016

@author: dongyu
"""

from netCDF4 import Dataset, num2date
import numpy as np
import matplotlib.pyplot as plt
import utm
from datetime import datetime
from scipy import interpolate
from mpl_toolkits.basemap import Basemap
from interpXYZ import interpXYZ
import pdb



class UV_interpolation(object):
    """
    This class is used to interpolate SUNTANS and ROMS velocity field onto
    the blended grid
    """
    def __init__(self,grid_file, sun_output, roms_file, blended_file, **kwargs):
        self.__dict__.update(kwargs)
        
        #### Read grid data and velocity data ####
        #self.readFile('blended_grid.nc','GalvCoarse6.nc','hiroms_ss_rho.nc')
        self.readFile(grid_file, sun_output, roms_file)
        
        #### Temporally and spatially interpolate SUNTANS velocity on blended grid ####
        #self.interp_sun()
        #self.interp_roms()
        #self.combine()
        #self.smooth_interp()
        #self.write_file('blendeduv.nc')


        self.write_file(blended_file)
        
        
        
    def readFile(self, blended_grid, suntans_output, roms_output):
        """
        functions used to read grid file and suntans output file
        """
        #### step 1) SUNTANS output
        #filename='GalvCoarse_0000.nc'
        nc = Dataset(suntans_output,'r')
        print "#### Reading SUNTANS output file !!!!     ####\n"
        #print nc

        self.uc=nc.variables['uc'][:][:,0,:]
        self.vc=nc.variables['vc'][:][:,0,:]
        xi=nc.variables['xv'][:]
        yi=nc.variables['yv'][:]
        timei=nc.variables['time']
        self.stime = num2date(timei[:],timei.units)   ## SUNTANS time: seconds since 1990-01-01
        
        
        #### step 2) ROMS output
        nc2 = Dataset(roms_output,'r')
        print "#### Reading ROMS output file !!!!        ####\n"
        lon0=nc2.variables['lon'][:]
        lat0=nc2.variables['lat'][:]
        self.mask0=nc2.variables['mask'][:]
        self.u_roms=nc2.variables['water_u'][:]
        self.v_roms=nc2.variables['water_v'][:]
        ftime=nc2.variables['time']
        rtime = num2date(ftime[:],ftime.units)  ## ROMS time: seconds since 1970-01-01        
        (y0,x0)=lon0.shape
        
        #### subset roms time, note that SUNTANS time period is shorted than ROMS time period
        t0 = self.stime[0]
        t1 = self.stime[-1]
        self.ind0 = self.findNearest(t0,rtime)
        self.ind1 = self.findNearest(t1,rtime)
        self.time_ss = rtime[self.ind0:self.ind1+1] 
        
        xroms0 = np.zeros_like(lon0)
        yroms0 = np.zeros_like(lat0)
        for i in range(y0):
            for j in range(x0):
                (yroms0[i,j],xroms0[i,j])=utm.from_latlon(lat0[i,j],lon0[i,j])[0:2]
        
        #### step 3) Blended grid
        #filename='blended_grid.nc'
        nc1 = Dataset(blended_grid,'r')
        print "#### Reading curvilinear blended grid !!!!####\n"
        #print nc1
        
        xr = nc1.variables['xr'][:]
        yr = nc1.variables['yr'][:]
        lon = nc1.variables['lon_rho'][:]
        lat = nc1.variables['lat_rho'][:]
        mask = nc1.variables['mask_rho'][:]
        
        xroms = np.zeros_like(lon)
        yroms = np.zeros_like(lat)
        (y,x) = lon.shape
        for i in range(y):
            for j in range(x):
                (yroms[i,j],xroms[i,j])=utm.from_latlon(lat[i,j],lon[i,j])[0:2]
        
        #### subset ROMS grid for interpolation ####
        def findNearset(x,y,lon,lat):
            """
            Return the J,I indices of the nearst grid cell to x,y
            """
                        
            dist = np.sqrt( (lon - x)**2 + (lat - y)**2)
            
            return np.argwhere(dist==dist.min())
                
        SW=utm.to_latlon(xi.min(),yi.min(),15,'R')  ###(lat, lon)
        NE=utm.to_latlon(xi.max(),yi.max(),15,'R')  
        
        #SW=utm.from_latlon(27.95,-95.2)[0:2]
        #NE=utm.from_latlon(30.0, -94.25)[0:2]
        
        #### step 4) searching for the index of the subset domain for interpolation
        #### Hard coded to narrow this domain
        ind = findNearset(SW[1], SW[0], lon, lat)
        J0=ind[0][0] + 15
        I0=ind[0][1] + 15
        
        ind = findNearset(NE[1], NE[0], lon, lat)
        J1=ind[0][0] + 25
        I1=ind[0][1] - 35
        
        yss = yroms[J0:J1,I0:I1]  ##subset x,y
        xss = xroms[J0:J1,I0:I1]
        self.maskss = mask[J0:J1,I0:I1]
        
        #pdb.set_trace()
        #### Step 5) Prepare the grid variables for the SUNTANS interpolation class
        xy_sun = np.vstack((yi.ravel(),xi.ravel())).T   ## SUNTANS grid, xi: latitude, yi: longitude        
        xy_new = np.vstack((xss[self.maskss==1],yss[self.maskss==1])).T   ## blended grid
        
        self.Fuv = interpXYZ(xy_sun,xy_new, method='idw')
        
        #### define spatial and length scales ####
        self.Nt = self.stime.shape[0]
        (self.X, self.Y) = xss.shape
        
        
        #### step 6) Prepare for interpolate original ROMS velocity
        SW0 = (lat0.min(),lon0.min())
        NE0 = (lat0.max(),lon0.max())
        ind0 = findNearset(SW0[1], SW0[0], lon, lat)
        JJ0=ind0[0][0] -40
        II0=ind0[0][1] 
        
        ind0 = findNearset(NE0[1], NE0[0], lon, lat)
        JJ1=ind0[0][0] +25
        II1=ind0[0][1] 
        
        yss0 = yroms[JJ0:JJ1,II0:II1]  ##subset x,y for ROMS velocity
        xss0 = xroms[JJ0:JJ1,II0:II1]
        self.maskss0 = mask[JJ0:JJ1,II0:II1]
        
        #### step 7) Prepare the grid variables for the SUNTANS interpolation class
        xy_roms = np.vstack((xroms0[self.mask0==1],yroms0[self.mask0==1])).T
        xy_new0 = np.vstack((xss0[self.maskss0==1],yss0[self.maskss0==1])).T   ## blended grid
        
        self.Fuv0 = interpXYZ(xy_roms,xy_new0, method='idw')
        
        #### define spatial and length scales ####
        self.Nt0 = rtime.shape[0]
        (self.X0, self.Y0) = xss0.shape
        self.lon0=lon[JJ0:JJ1,II0:II1]
        self.lat0=lat[JJ0:JJ1,II0:II1]
        
        #### step 8) define the index of SUNTANS in sub-domain of ROMS
        self.JJJ0 = J0-JJ0-1
        self.JJJ1 = J1-JJ0
        
        self.III0 = I0-II0
        self.III1 = I1-II0
        
        #### the new time is the ROMS output time ####
        self.time = ftime[:]
        
        
        #pdb.set_trace()
        
        
        
    def interp_sun(self):
        """
        interpolate SUNTANS velocity onto blended grid spatially and temporally 
        Performs the interpolation in this order:
            1) Interpolate onto horizontal coordinates
            2) Interpolate onto time coordinates
        """
        #### Initialize the output array @ subset roms grid ####
        uout=np.zeros([self.Nt, self.X, self.Y])
        vout=np.zeros([self.Nt, self.X, self.Y])            
        
        print 'Spatially interpolating SUNTANS u, v'
        #### Loop through each time step ####
        for tstep in range(self.Nt):
            utem = self.Fuv(self.uc[tstep,:])
            vtem = self.Fuv(self.vc[tstep,:])
            uout[tstep][self.maskss==1] = utem
            vout[tstep][self.maskss==1] = vtem
        
        # interpolate temporally
        print 'Temperally interpolating SUNTANS u, v'
        troms = self.SecondsSince(self.time_ss)
        tsuntans = self.SecondsSince(self.stime)
        
        Ft = interpolate.interp1d(tsuntans,uout,axis=0,kind='linear',bounds_error=False)
        uss = Ft(troms)    ## temporally and spatially interpolated subset u velocity
        
        Ft = interpolate.interp1d(tsuntans,vout,axis=0,kind='linear',bounds_error=False)
        vss = Ft(troms)    ## temporally and spatially interpolated subset v velocity        
        
        plt.figure(1)
        plt.pcolor(uout[20,:,:])
        plt.ylim((0,self.X))
        plt.xlim((0,self.Y))
        #plt.show()
        
        return uss, vss
          
        
    def interp_roms(self):
        """
        interpolate ROMS velocity onto blended grid spatially and temporally
        """
        #### Initialize the output array @ subset roms grid ####
        uout=np.zeros([self.Nt0, self.X0, self.Y0])
        vout=np.zeros([self.Nt0, self.X0, self.Y0])            
        
        print 'Spatially interpolating ROMS u, v'
        #### Loop through each time step ####
        for tstep in range(self.Nt0):
            utem = self.Fuv0(self.u_roms[tstep,:,:][self.mask0==1].flatten())
            vtem = self.Fuv0(self.v_roms[tstep,:,:][self.mask0==1].flatten())
            #pdb.set_trace()
            uout[tstep,:,:][self.maskss0==1] = utem
            vout[tstep,:,:][self.maskss0==1] = vtem
                 
        
        basemap = Basemap(projection='merc',llcrnrlat=self.lat0.min(),urcrnrlat=self.lat0.max(), \
                    llcrnrlon=self.lon0.min(),urcrnrlon=self.lon0.max(),resolution='i')
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        
        basemap.drawcoastlines()
        basemap.fillcontinents()
        basemap.drawcountries()
        basemap.drawstates()
        x_rho, y_rho = basemap(self.lon0, self.lat0)

        basemap.pcolormesh(x_rho, y_rho, uout[-2,:,:],vmin=-0.55588271196095485,vmax=0.94814157426936685)        
        #plt.show()

        return uout, vout        
        
        
    def combine(self):
        """
        insert suntans velocity into roms
        """
        usun, vsun = self.interp_sun()      
        uroms, vroms = self.interp_roms()
        uroms[self.ind0:self.ind1+1,self.JJJ0:self.JJJ1,self.III0:self.III1]=usun
        vroms[self.ind0:self.ind1+1,self.JJJ0:self.JJJ1,self.III0:self.III1]=vsun
        
        basemap = Basemap(projection='merc',llcrnrlat=self.lat0.min(),urcrnrlat=self.lat0.max(), \
                    llcrnrlon=self.lon0.min(),urcrnrlon=self.lon0.max(),resolution='i')
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        
        basemap.drawcoastlines()
        basemap.fillcontinents()
        basemap.drawcountries()
        basemap.drawstates()
        x_rho, y_rho = basemap(self.lon0, self.lat0)

        basemap.pcolormesh(x_rho, y_rho, uroms[-2,:,:], vmin=uroms.min(),vmax=uroms.max())        
        #plt.show()
        
        return uroms, vroms
        
    def smooth_interp(self):
        """
        re-interpolate the new ROMS velocity field to make the boundary variation smooth
        """
        #### 1) Prepare for ROMS grid
        xroms=np.zeros_like(self.lon0)  ## xroms: longitude, yroms: latitude
        yroms=np.zeros_like(self.lat0)
        (y,x)=self.lon0.shape
        for i in range(y):
            for j in range(x):
                (yroms[i][j],xroms[i][j])=utm.from_latlon(self.lat0[i][j],self.lon0[i][j])[0:2]
        
        xy_roms = np.vstack((xroms[self.maskss0==1],yroms[self.maskss0==1])).T
        Fuv = interpXYZ(xy_roms,xy_roms, method='kriging')
        
        uroms, vroms = self.combine()
        for tstep in range(self.time_ss.shape[0]):
            utem=np.zeros_like(xroms)
            utem[self.maskss0==1]=Fuv(uroms[self.ind0+tstep,:,:][self.maskss0==1])
            uroms[self.ind0+tstep,:,:]=utem
            
            vtem=np.zeros_like(xroms)
            vtem[self.maskss0==1]=Fuv(vroms[self.ind0+tstep,:,:][self.maskss0==1])
            vroms[self.ind0+tstep,:,:]=vtem
        
        basemap = Basemap(projection='merc',llcrnrlat=self.lat0.min(),urcrnrlat=self.lat0.max(), \
                    llcrnrlon=self.lon0.min(),urcrnrlon=self.lon0.max(),resolution='i')
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        
        basemap.drawcoastlines()
        basemap.fillcontinents()
        basemap.drawcountries()
        basemap.drawstates()
        x_rho, y_rho = basemap(self.lon0, self.lat0)

        basemap.pcolormesh(x_rho, y_rho, uroms[-2,:,:], vmin=uroms.min(),vmax=uroms.max())        
        plt.show()        
        
        #pdb.set_trace()
    
    def write_file(self,outfile):
        """
        write the netcdf file using the updated velocity field
        """
        print "dumping the interpolating velocity fields into %s!!!\n"%outfile
        uroms, vroms = self.combine()

        y,x=uroms.shape[1:]
        nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        nc.grid_type = 'curvilinear'
        nc.set_fill_off()
        ####Create dimensions####
        nc.createDimension('x', x) 
        nc.createDimension('y', y)
        nc.createDimension('xc', x) 
        nc.createDimension('yc', y)
        nc.createDimension('time', None)   ##unlimited dimension
        nc.close()
        
        ####adding variables####
        
        self.create_nc_var(outfile,'time',('time'),{'units':'seconds since 1970-01-01 00:00:00'})
        self.create_nc_var(outfile,'lon',('y','x'),{'long_name':'Longitude',\
            'units':'degrees_east','standard_name':'longitude'})
        self.create_nc_var(outfile,'lat',('y','x'),{'long_name':'Latitude',\
            'units':'degrees_north','standard_name':'latitude'})
        self.create_nc_var(outfile,'water_u',('time','yc','xc'),{ \
            'units':'meter second-1'})
        self.create_nc_var(outfile,'water_v',('time','yc','xc'),{ \
            'units':'meter second-1'})   
        self.create_nc_var(outfile,'mask',('yc','xc'),{})
        
        
        ######Now writting the variables######
        nc = Dataset(outfile,'a')
        nc.variables['time'][:] = self.time
        nc.variables['lon'][:] = self.lon0
        nc.variables['lat'][:] = self.lat0
        nc.variables['water_u'][:] = uroms
        nc.variables['water_v'][:] = vroms
        nc.variables['mask'] = self.maskss0
        print "Ending writing variables into %s netcdf file !!!"%outfile
        nc.close()
        
        
        
    def findNearest(self,t,timevec):
        """
        Return the index from timevec the nearest time point to time, t. 
        
        """
        tnow = self.SecondsSince(t)
        tvec = self.SecondsSince(timevec)
        
        #tdist = np.abs(tnow - tvec)
        
        #idx = np.argwhere(tdist == tdist.min())
        
        #return int(idx[0])
        return np.searchsorted(tvec,tnow)[0]
        
        
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
        
    def create_nc_var(self,outfile, name, dimensions, attdict, dtype='f8',zlib=False,complevel=0):
        
        nc = Dataset(outfile, 'a')
        tmp=nc.createVariable(name, dtype, dimensions,zlib=zlib,complevel=complevel)
        for aa in attdict.keys():
            tmp.setncattr(aa,attdict[aa])
        #nc.variables[name][:] = var	
        nc.close()

########################### For testing #############################
if __name__ == "__main__":
    UV_interpolation()
    print "end"
