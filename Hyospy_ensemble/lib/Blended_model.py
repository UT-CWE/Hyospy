# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 10:41:02 2016

@author: dongyu
"""

import numpy as np
import utm
import string
import os
import glob
from datetime import datetime
from matplotlib.collections import PolyCollection
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from interpXYZ import interpXYZ
from romsio import roms_subset
from SUNTANSio import SUNTANS_download
import pdb

class blend(object):
    """
    download the data from the server and blend model result
    superposition
    """

    def __init__(self, starttime, endtime, **kwargs):
        self.__dict__.update(kwargs)
        
        timelims = (starttime.replace("-", "")+'0000', endtime.replace("-", "")+'0000')

        #### SUNTANS file ####
        basedir = os.getcwd()
        sun_dir = basedir + '/SUNTANS/rundata'
        os.chdir(sun_dir)
        ncfiles = []
        for ff in glob.glob("GalvCoarse_0*"):
            ncfiles.append(ff)
	
        os.chdir(basedir)

        nc_suntans = []
        for f in ncfiles:
            nc_suntans.append('%s/%s'%(sun_dir, f))	    
        
        self.sun = SUNTANS_download(nc_suntans, timelims)
	
        #### ROMS file ####
        grdfile = 'DATA/txla_grd_v4_new.nc'
        #nc_roms = ['http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg']
        nc_roms = ['DATA/txla_subset_HIS.nc']

        bbox = [-98.96,-87.80,26.00,32.24]  ## the whole domain
        self.roms = roms_subset(nc_roms,bbox,timelims,gridfile=grdfile)
     
        ## Read weight variables ##
        nc_w = Dataset('DATA/weights.nc', 'r')
        self.w_sun = nc_w.variables['w_sun'][:]
        self.w_roms = nc_w.variables['w_roms'][:]
        nc_w.close()
          
        #self.model_velocity()
              

    def model_velocity(self):
        """
        Read SUNTANS and ROMS velocity
        """
        
        #### ROMS data ####
        self.mask_roms = self.roms.mask_rho[0:-1,0:-1]
        self.lon_rho = self.roms.lon_rho[0:-1,0:-1]
        self.lat_rho = self.roms.lat_rho[0:-1,0:-1]
        ang = self.roms.angle[0:-1,0:-1]
        
        #### load blended grid ####
        self.blended_grid()

        #### create blended_uv.nc file to dump the data ####
        self.Writefile('DATA/blended_uv.nc')        
        
        #### empty matrix for SUNTANS and ROMS
        self.usun = np.zeros([self.sun.Nt, self.xss.shape[0], self.xss.shape[1]])
        self.vsun = np.zeros([self.sun.Nt, self.xss.shape[0], self.xss.shape[1]])
        self.uroms = np.zeros([self.roms.Nt, self.xss0.shape[0], self.yss0.shape[1]])
        self.vroms = np.zeros([self.roms.Nt, self.xss0.shape[0], self.yss0.shape[1]])
        
        self.eta_sun = np.zeros([self.sun.Nt, self.xss.shape[0], self.xss.shape[1]])
        self.eta_roms = np.zeros([self.roms.Nt, self.xss0.shape[0], self.yss0.shape[1]])
        pdb.set_trace()
        #### Loop through every time step ####
        if self.sun.Nt == self.roms.Nt:            
            for ii in range(0, self.roms.Nt):
                ## SUNTANS velocity
                self.sun.ReadData(ii)
                uc = self.sun.uc[2,:]
                vc = self.sun.vc[2,:]
                eta = self.sun.eta[:]
                
                ## plot suntans velocity at time ii for testing                
                #self.plot_suntans(ii, eta, self.sun.timei[ii])
                ## interp SUNTANS
                print 'Interpolating SUNTANS velocity at time: %s...'%datetime.strftime(self.sun.timei[ii],'%Y-%m-%d %H:%M:%S')
                
                self.usun[ii,:,:][(self.maskss==2)|(self.maskss==3)|(self.maskss==4)|(self.maskss==5)|(self.maskss==6)],  \
                self.vsun[ii,:,:][(self.maskss==2)|(self.maskss==3)|(self.maskss==4)|(self.maskss==5)|(self.maskss==6)],  \
                self.eta_sun[ii,:,:][(self.maskss==2)|(self.maskss==3)|(self.maskss==4)|(self.maskss==5)|(self.maskss==6)], \
                = self.interp_sun(uc, vc, eta)
                
                ## ROMS velocity
                self.roms.ReadData(ii)
                u_roms = self.shrink(self.roms.u[-1,:,:], self.mask_roms.shape)
                v_roms = self.shrink(self.roms.v[-1,:,:], self.mask_roms.shape)                
                u_roms, v_roms =  self.rot2d(u_roms, v_roms, ang)
                
                eta_roms = self.shrink(self.roms.zeta[:,:], self.mask_roms.shape)   
                
                ## plot roms velocity at time ii for testing
                #self.plot_roms(ii, eta_roms, self.roms.time[ii])
                ## interp ROMS
                print 'Interpolating ROMS velocity at time: %s...'%datetime.strftime(self.roms.time[ii],'%Y-%m-%d %H:%M:%S')
                
                
                self.uroms[ii,:,:], self.vroms[ii,:,:], self.eta_roms[ii,:,:] \
                                = self.interp_roms_new(u_roms, v_roms, eta_roms)
                

                ## insert SUNTANS to blended grid
                self.uroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==2] \
                = self.usun[ii,:,:][self.maskss==2]
                self.vroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==2] \
                = self.vsun[ii,:,:][self.maskss==2]
                self.eta_roms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==2] \
                = self.eta_sun[ii,:,:][self.maskss==2]
                
                ## buffer domain
                print 'Buffer domain smoothing at time: %s...'%datetime.strftime(self.roms.time[ii],'%Y-%m-%d %H:%M:%S')
                #### u-velocity ####
                self.uroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==6] \
                    = self.usun[ii,:,:][self.maskss==6] * self.w_sun[self.maskss0==6] + \
                      self.uroms[ii,:,:][self.maskss0==6] * self.w_roms[self.maskss0==6]
                self.uroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==3] \
                    = self.usun[ii,:,:][self.maskss==3] * self.w_sun[self.maskss0==3] + \
                      self.uroms[ii,:,:][self.maskss0==3] * self.w_roms[self.maskss0==3]
                self.uroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==4] \
                    = self.usun[ii,:,:][self.maskss==4] * self.w_sun[self.maskss0==4] + \
                      self.uroms[ii,:,:][self.maskss0==4] * self.w_roms[self.maskss0==4]
                # superposition      
                self.uroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==5] \
                    = self.usun[ii,:,:][self.maskss==5] * self.w_sun[self.maskss0==5] + \
                      self.uroms[ii,:,:][self.maskss0==5] * self.w_roms[self.maskss0==5]
                      
                #### v-velocity ####
                self.vroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==6] \
                    = self.vsun[ii,:,:][self.maskss==6] * self.w_sun[self.maskss0==6] + \
                      self.vroms[ii,:,:][self.maskss0==6] * self.w_roms[self.maskss0==6]
                self.vroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==3] \
                    = self.vsun[ii,:,:][self.maskss==3] * self.w_sun[self.maskss0==3] + \
                      self.vroms[ii,:,:][self.maskss0==3] * self.w_roms[self.maskss0==3]
                self.vroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==4] \
                    = self.vsun[ii,:,:][self.maskss==4] * self.w_sun[self.maskss0==4] + \
                      self.vroms[ii,:,:][self.maskss0==4] * self.w_roms[self.maskss0==4]
                # superposition
                self.vroms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==5] \
                    = self.vsun[ii,:,:][self.maskss==5] * self.w_sun[self.maskss0==5] + \
                      self.vroms[ii,:,:][self.maskss0==5] * self.w_roms[self.maskss0==5]
                      
                #### eta ####
                self.eta_roms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==6] \
                    = self.eta_sun[ii,:,:][self.maskss==6] * self.w_sun[self.maskss0==6] + \
                      self.eta_roms[ii,:,:][self.maskss0==6] * self.w_roms[self.maskss0==6]
                self.eta_roms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==3] \
                    = self.eta_sun[ii,:,:][self.maskss==3] * self.w_sun[self.maskss0==3] + \
                      self.eta_roms[ii,:,:][self.maskss0==3] * self.w_roms[self.maskss0==3]
                self.eta_roms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==4] \
                    = self.eta_sun[ii,:,:][self.maskss==4] * self.w_sun[self.maskss0==4] + \
                      self.eta_roms[ii,:,:][self.maskss0==4] * self.w_roms[self.maskss0==4]
                      
                self.eta_roms[ii,self.JJJ0:self.JJJ1,self.III0:self.III1][self.maskss==5] \
                    = self.eta_sun[ii,:,:][self.maskss==5] * self.w_sun[self.maskss0==5] + \
                      self.eta_roms[ii,:,:][self.maskss0==5] * self.w_roms[self.maskss0==5]
                                            
                #### Test Plotting ####
#                basemap = Basemap(projection='merc',llcrnrlat=self.lat0.min(),urcrnrlat=self.lat0.max(), \
#                    llcrnrlon=self.lon0.min(),urcrnrlon=self.lon0.max(),resolution='i')
#                fig1 = plt.figure()
#                ax = fig1.add_subplot(111)
#                
#                basemap.drawcoastlines()
#                basemap.fillcontinents()
#                basemap.drawcountries()
#                basemap.drawstates()
#                x_rho, y_rho = basemap(self.lon0, self.lat0)
#                
#                basemap.pcolor(x_rho, y_rho, self.uroms[ii,:,:], vmin=-0.7,vmax=0.7) 
#                plt.title('Blended velocity at time: %s...'%datetime.strftime(self.roms.time[ii],'%Y-%m-%d %H:%M:%S'))
#                #plt.savefig(os.getcwd()+'/figures/blended_u'+str(ii)+'.png')
#                plt.show()
                
                self.Writedata(ii, self.uroms[ii,:,:], self.vroms[ii,:,:], self.eta_roms[ii,:,:])
                
                #pdb.set_trace()
                
        else:
            print "Time length of SUNTANS output is different from that of ROMS output!!!"
            
        print '##################\nDone!\n##################'
    
    
    def interp_sun(self, usuntans, vsuntans, eta_sun):
        """
        function that interpolate SUNTANS velocity onto blended grid
        """
        ## 1) Prepare the interplation class
        #### SUNTANS
        xy_sun = np.vstack((self.sun.yv.ravel(),self.sun.xv.ravel())).T
        xy_new = np.vstack((self.xss[(self.maskss==2)|(self.maskss==3)|(self.maskss==4)|(self.maskss==5)|(self.maskss==6)], \
                                self.yss[(self.maskss==2)|(self.maskss==3)|(self.maskss==4)|(self.maskss==5)|(self.maskss==6)])).T    # blended grid        
        Fuv = interpXYZ(xy_sun, xy_new, method='idw')
        
        return Fuv(usuntans), Fuv(vsuntans), Fuv(eta_sun)
                
    def interp_roms(self, uroms, vroms):
        """
        function that interpolate ROMS velocity onto blended grid
        """
        #### ROMS 
        xroms0 = np.zeros_like(self.lon_rho)
        yroms0 = np.zeros_like(self.lat_rho)
        for i in range(self.lon_rho.shape[0]):
            for j in range(self.lat_rho.shape[1]):
                (yroms0[i,j],xroms0[i,j])=utm.from_latlon(self.lat_rho[i,j],self.lon_rho[i,j])[0:2]
                
        xy_roms = np.vstack((xroms0[self.mask_roms==1],yroms0[self.mask_roms==1])).T
        xy_new0 = np.vstack((self.xss0[(self.maskss0==1)|(self.maskss0==3)|(self.maskss0==4)], \
                                self.yss0[(self.maskss0==1)|(self.maskss0==3)|(self.maskss0==4)])).T
                                
        Fuv0 = interpXYZ(xy_roms,xy_new0, method='idw')
        
        return Fuv0(uroms[:,:][(self.mask_roms==1)].flatten()), Fuv0(vroms[:,:][(self.mask_roms==1)].flatten())

    def interp_roms_new(self, uroms, vroms, eta_roms):
        """
        Step 1: divide roms velocity field into two parts
        Step 2: interpolate each
        Step 3: combine them together
        """
        #### subset Blend_grid for interpolation ####
        def findNearset(x,y,lon,lat):
            """
            Return the J,I indices of the nearst grid cell to x,y
            """
                        
            dist = np.sqrt( (lon - x)**2 + (lat - y)**2)
            
            return np.argwhere(dist==dist.min())        
        
        xroms0 = np.zeros_like(self.lon_rho)
        yroms0 = np.zeros_like(self.lat_rho)
        for i in range(self.lon_rho.shape[0]):
            for j in range(self.lat_rho.shape[1]):
                (yroms0[i,j],xroms0[i,j])=utm.from_latlon(self.lat_rho[i,j],self.lon_rho[i,j])[0:2]
                
        xy_roms = np.vstack((xroms0[self.mask_roms==1],yroms0[self.mask_roms==1])).T
        J0 = 50; J1 = -1; I0 = 142; I1 = 580
        
        SW = [self.lat0[J0,I0], self.lon0[J0,I0]]
        NE = [self.lat0[J1,I1], self.lon0[J1,I1]]
        
        ind = findNearset(SW[1], SW[0], self.lon_rho, self.lat_rho)
        JJ0=ind[0][0] 
        II0=ind[0][1] 
        
        ind = findNearset(NE[1], NE[0], self.lon_rho, self.lat_rho)
        JJ1=ind[0][0]
        II1=ind[0][1]
        
        mask_old = self.mask_roms.copy()
        mask_old[JJ0:JJ1, II0:II1][mask_old[JJ0:JJ1,II0:II1]==1] = 2
        
        xy_roms1 = np.vstack((xroms0[mask_old==2], \
                                yroms0[mask_old==2])).T
        xy_roms2 = np.vstack((xroms0[mask_old==1], \
                                yroms0[mask_old==1])).T
        

        mask = self.maskss0.copy()
        mask[(mask==2)|(mask==3)|(mask==4)|(mask==5)|(mask==6)] = 1        
        
        ## mask1==2 is in the core area; mask1==1 is in other broader area
        mask1 = mask.copy()
        mask1[J0:, I0:I1][mask1[J0:,I0:I1]==1] = 2
        
                
        xy_new1 = np.vstack((self.xss0[mask1==2], self.yss0[mask1==2])).T
        xy_new2 = np.vstack((self.xss0[mask1==1], self.yss0[mask1==1])).T

        uout = np.zeros([self.xss0.shape[0], self.xss0.shape[1]])
        vout = np.zeros([self.xss0.shape[0], self.xss0.shape[1]])
        eta_out = np.zeros([self.xss0.shape[0], self.xss0.shape[1]])
        
        
        Fuv1 = interpXYZ(xy_roms1, xy_new1, method='idw')
        uout[:,:][mask1==2], vout[:,:][mask1==2], eta_out[:,:][mask1==2] = \
            Fuv1(uroms[:,:][(mask_old==2)].flatten()), Fuv1(vroms[:,:][(mask_old==2)].flatten()), \
                Fuv1(eta_roms[:,:][(mask_old==2)].flatten())
            
        Fuv2 = interpXYZ(xy_roms2, xy_new2, method='idw')
        uout[:,:][mask1==1], vout[:,:][mask1==1], eta_out[:,:][mask1==1] = \
            Fuv2(uroms[:,:][(mask_old==1)].flatten()), Fuv2(vroms[:,:][(mask_old==1)].flatten()), \
                Fuv2(eta_roms[:,:][(mask_old==1)].flatten())
	
        return uout, vout, eta_out
        
        
    def blended_grid(self):
        """
        deal with the blended grid
        """
        
        nc = Dataset('DATA/blended_grid_new.nc', 'r') 
        #print nc
        
        lon = nc.variables['lon_rho'][:]
        lat = nc.variables['lat_rho'][:]
        mask = nc.variables['mask_rho'][:]
        x_bg = np.zeros_like(lon)
        y_bg = np.zeros_like(lat)
        for i in range(lon.shape[0]):
            for j in range(lon.shape[1]):
                (y_bg[i,j],x_bg[i,j])=utm.from_latlon(lat[i,j],lon[i,j])[0:2]
        
        #### subset Blend_grid for interpolation ####
        def findNearset(x,y,lon,lat):
            """
            Return the J,I indices of the nearst grid cell to x,y
            """
                        
            dist = np.sqrt( (lon - x)**2 + (lat - y)**2)
            
            return np.argwhere(dist==dist.min())

        #### Step 1) subset for SUNTANS interpolation
        NE = (29.868007, -94.175217) 
        SW = (28.361303, -95.073081) 
        
        #### searching for the index of the subset domain for interpolation
        ind = findNearset(SW[1], SW[0], lon, lat)
        J0=ind[0][0] 
        I0=ind[0][1] 
        
        ind = findNearset(NE[1], NE[0], lon, lat)
        J1=ind[0][0] +22 #-2
        I1=ind[0][1]
        
        self.yss = y_bg[J0:J1,I0:I1]  ##subset x,y
        self.xss = x_bg[J0:J1,I0:I1]
        self.maskss = mask[J0:J1,I0:I1]    
        
        #### Step 2) subset for ROMS interpolation
        #### Prepare for interpolate original ROMS velocity
        SW0 = (self.lat_rho[0,0], self.lon_rho[0,0])
        NE0 = (self.lat_rho[-1,-1], self.lon_rho[-1,-1])
        ind0 = findNearset(SW0[1], SW0[0], lon, lat)
        JJ0=ind0[0][0] 
        II0=ind0[0][1] 
        
        ind0 = findNearset(NE0[1], NE0[0], lon, lat)
        JJ1=ind0[0][0] + 6
        II1=ind0[0][1]
        
        #self.yss0 = y_bg[JJ0:JJ1,II0:II1]  ##subset x,y for ROMS velocity
        #self.xss0 = x_bg[JJ0:JJ1,II0:II1]
        #self.maskss0 = mask[JJ0:JJ1,II0:II1]
        #self.lon0 = lon[JJ0:JJ1,II0:II1]
        #self.lat0 = lat[JJ0:JJ1,II0:II1]
        
        self.yss0 = y_bg  ##subset x,y for ROMS velocity
        self.xss0 = x_bg
        self.maskss0 = mask
        self.lon0 = lon
        self.lat0 = lat
        
        self.lon_rho0 = nc.variables['lon_rho'][:]
        self.lat_rho0 = nc.variables['lat_rho'][:]
        self.angle_rho0 = nc.variables['angle_rho'][:]
        self.pm0 = nc.variables['pm'][:]
        self.pn0 = nc.variables['pn'][:]
        self.lon_vert0 = nc.variables['lon_vert'][:]
        self.lat_vert0 = nc.variables['lat_vert'][:]
        self.h0 = nc.variables['h'][:]

        self.lat_u = nc.variables['lat_u'][:]
        self.lon_u = nc.variables['lon_u'][:]
        self.mask_u = nc.variables['mask_u'][:]
        self.lat_v = nc.variables['lat_v'][:]
        self.lon_v = nc.variables['lon_v'][:]
        self.mask_v = nc.variables['mask_v'][:]
        self.angle = nc.variables['angle'][:]
        
        ## subset weight variables ##
        #self.w_sun = self.w_sun[JJ0:JJ1,II0:II1]
        #self.w_roms = self.w_roms[JJ0:JJ1,II0:II1]
        
        #self.JJJ0 = J0-JJ0
        #self.JJJ1 = J1-JJ0
        
        #self.III0 = I0-II0
        #self.III1 = I1-II0

        self.JJJ0 = J0
        self.JJJ1 = J1
        
        self.III0 = I0
        self.III1 = I1        
        #pdb.set_trace()

        
          
    def Writefile(self, outfile, verbose=True):
        """
        Write the new interpolated velocity file into a netcdf file
        """
        self.outfile = outfile
        
        # Write the grid to file
        nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        nc.Description = 'blended velocity file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()
        nc.type = 'blended HIS file'
        #pdb.set_trace()
        nc.createDimension('xr', self.lon0.shape[1])
        nc.createDimension('yr', self.lon0.shape[0])    
        nc.createDimension('ocean_time', None)
        
        nc.createDimension('xpsi', self.lon0.shape[1]-1)
        nc.createDimension('ypsi', self.lon0.shape[0]-1) 

        nc.createDimension('xvert', self.lon_vert0.shape[1])
        nc.createDimension('yvert', self.lon_vert0.shape[0])
        nc.createDimension('zl', 1)
        nc.createDimension('zlp1', 2)
        
        def write_nc_var(var, name, dimensions, units=None):
            nc.createVariable(name, 'f8', dimensions)
            if units is not None:
                nc.variables[name].units = units
            nc.variables[name][:] = var
            if verbose:
                print ' ... wrote ', name
                
        def create_nc_var(name, dimensions, units=None):
            nc.createVariable(name, 'f8', dimensions)
            if units is not None:
                nc.variables[name].units = units
            if verbose:
                print ' ... wrote ', name
                
        # Grid variables
        write_nc_var(self.lon0, 'lon', ('yr','xr'), 'meters')
        write_nc_var(self.lat0, 'lat', ('yr','xr'), 'meters')
        mask = self.maskss0.copy()
        mask[(mask==2)|(mask==3)|(mask==4)] = 1
        write_nc_var(mask, 'mask', ('yr','xr'))

        write_nc_var(mask, 'mask_rho', ('yr','xr'))
        write_nc_var(self.lon_rho0, 'lon_rho', ('yr','xr'), 'meters')
        write_nc_var(self.lat_rho0, 'lat_rho', ('yr','xr'), 'meters')
        write_nc_var(self.lon_u, 'lon_u', ('yr','xpsi'), 'meters')
        write_nc_var(self.lat_u, 'lat_u', ('yr','xpsi'), 'meters')
        write_nc_var(self.mask_u, 'mask_u', ('yr','xpsi'))
        write_nc_var(self.lon_v, 'lon_v', ('ypsi','xr'), 'meters')
        write_nc_var(self.lat_v, 'lat_v', ('ypsi','xr'), 'meters')
        write_nc_var(self.mask_v, 'mask_v', ('ypsi','xr'))
        write_nc_var(self.angle, 'angle', ('yvert','xvert'))
        
        write_nc_var(self.angle_rho0, 'angle_rho', ('yr','xr'))
        write_nc_var(self.lon_vert0, 'lon_vert', ('yvert','xvert'), 'meters')
        write_nc_var(self.lat_vert0, 'lat_vert', ('yvert','xvert'), 'meters')
        write_nc_var(self.h0, 'h', ('yr','xr'))
        write_nc_var(self.pm0, 'pm', ('yr','xr'))
        write_nc_var(self.roms.hc, 'hc', ())
        write_nc_var(self.roms.Vstretching, 'Vstretching', ())
        write_nc_var(self.roms.Vtransform, 'Vtransform', ())
        write_nc_var(self.roms.theta_s, 'theta_s', ())
        write_nc_var(self.roms.theta_b, 'theta_b', ())
        write_nc_var(self.roms.s_rho[-1], 's_rho', ('zl'))
        write_nc_var(self.roms.Cs_r[-1], 'Cs_r', ('zl'))
        write_nc_var(self.roms.s_w[-3:-1], 's_w', ('zlp1'))
        write_nc_var(self.roms.Cs_w[-3:-1], 'Cs_w', ('zlp1'))
        
        # Create the data variables
        create_nc_var('ocean_time',('ocean_time'),'seconds since 1970-01-01 00:00:00')
        create_nc_var('u',('ocean_time','yr','xpsi'),'meter second-1')
        create_nc_var('v',('ocean_time','ypsi','xr'),'meter second-1')
        create_nc_var('zeta',('ocean_time','yr','xr'),'meter second-1')
        
        nc.close()
        
        
    def Writedata(self,tstep, u, v, eta):
        """
        funtion that writes the data at every time step into the variables
        """
        nc = Dataset(self.outfile, 'a')
        nc.variables['ocean_time'][tstep] = self.roms.ocean_time
        nc.variables['u'][tstep,:,:] = u[:,:-1]
        nc.variables['v'][tstep,:,:] = v[:-1,:]
        nc.variables['zeta'][tstep,:,:] = eta[:,:]
        
        nc.close()
   
    def plot_suntans(self, t, var, time):
        """
        function that plots SUNTANS velocity at one time step
        """
        timeformat = '%Y%m%d-%H%M'
        maxfaces = self.sun.cells.shape[1]
        basedir = os.getcwd() + '/grids'
        pointdata = self.readTXT(basedir+'/points.dat')
        xp1 = pointdata[:,0]
        yp1 = pointdata[:,1]
        #########################################################
        xp = np.zeros((self.sun.Nc,maxfaces+1))
        yp = np.zeros((self.sun.Nc,maxfaces+1))
            
        cells=self.sun.cells.copy()
        nfaces = 3*np.ones((self.sun.Nc,),np.int)        
        
        
        xp[:,:maxfaces]=xp1[cells]
        xp[range(self.sun.Nc),nfaces]=xp1[cells[:,0]]
        yp[:,:maxfaces]=yp1[cells]
        yp[range(self.sun.Nc),nfaces]=yp1[cells[:,0]]
        
        xy = np.zeros((maxfaces+1,2))
        def _closepoly(ii):
            nf=nfaces[ii]+1
            xy[:nf,0]=xp[ii,:nf]
            xy[:nf,1]=yp[ii,:nf]
            return xy[:nf,:].copy()

        cellxy= [_closepoly(ii) for ii in range(self.sun.Nc)]

        #clim=[w.min(),w.max()]
        clim=[-0.0010,0.0010]
        xlims=(self.sun.xv.min(),self.sun.xv.max())
        ylims=(self.sun.yv.min(),self.sun.yv.max())  
        
        fig = plt.figure(figsize=(10,8))
        axes = fig.add_subplot(111)
        collection = PolyCollection(cellxy,cmap='jet')
        collection.set_array(np.array(var[:]))
        collection.set_edgecolors('k')
        collection.set_linewidths(0.2)
        #collection.set_clim(vmin=clim[0],vmax=clim[1])
        collection.set_edgecolors(collection.to_rgba(np.array((var[:])))) 
        cbar = fig.colorbar(collection,orientation='vertical')
        axes.add_collection(collection)
        axes.set_xlim(xlims)
        axes.set_ylim(ylims)
        axes.set_aspect('equal')
        axes.set_xlabel('Easting [m]')
        axes.set_ylabel('Northing [m]')
        timestr = datetime.strftime(time, timeformat)
        plt.title('SUNTANS U velocity at %s'%timestr)
        plt.savefig(os.getcwd()+'/figures/suntans_eta_'+str(t)+'.png')
        #plt.show()
        
    def plot_roms(self, ii, var, time):
        """
        function that plots ROMS velocity at the time
        """
        
        basemap = Basemap(projection='merc',llcrnrlat=self.lat_rho.min(),urcrnrlat=self.lat_rho.max(), \
                    llcrnrlon=self.lon_rho.min(),urcrnrlon=self.lon_rho.max(),resolution='i')
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        
        basemap.drawcoastlines()
        basemap.fillcontinents()
        basemap.drawcountries()
        basemap.drawstates()
        x_rho, y_rho = basemap(self.lon_rho, self.lat_rho)
        
        basemap.pcolor(x_rho, y_rho, var[:,:], vmin=-0.7,vmax=0.7) 
        plt.title('ROMS U velocity at time: %s...'%datetime.strftime(time,'%Y-%m-%d %H:%M:%S'))
        plt.savefig(os.getcwd()+'/figures/roms_eta_'+str(ii)+'.png')
        #plt.show()
        
        
    def rot2d(self, x, y, ang):
        """
        rotate vectors by geometric angle
        This routine is part of Rob Hetland's OCTANT package:
        https://github.com/hetland/octant
        """
        xr = x*np.cos(ang) - y*np.sin(ang)
        yr = x*np.sin(ang) + y*np.cos(ang)
        return xr, yr
    
    
    def _int2str(self, num):
        """
        convert the month format
        input a integer;
        output a string;
        """
        if num<10:
            return '00%s'%str(num)
        elif 10<=num<100:
            return '0%s'%str(num)
        else:
            return '%s'%str(num)
            
    def shrink(self,a,b):
        """Return array shrunk to fit a specified shape by triming or averaging.
        
        a = shrink(array, shape)
        
        array is an numpy ndarray, and shape is a tuple (e.g., from
        array.shape). a is the input array shrunk such that its maximum
        dimensions are given by shape. If shape has more dimensions than
        array, the last dimensions of shape are fit.
        
        as, bs = shrink(a, b)
        
        If the second argument is also an array, both a and b are shrunk to
        the dimensions of each other. The input arrays must have the same
        number of dimensions, and the resulting arrays will have the same
        shape.
        
        This routine is part of Rob Hetland's OCTANT package:
            https://github.com/hetland/octant
            
        Example
        -------
        
        >>> shrink(rand(10, 10), (5, 9, 18)).shape
        (9, 10)
        >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
        [(5, 9, 10), (5, 9, 10)]   
        
        """

        if isinstance(b, np.ndarray):
            if not len(a.shape) == len(b.shape):
                raise Exception, \
                      'input arrays must have the same number of dimensions'
            a = self.shrink(a,b.shape)
            b = self.shrink(b,a.shape)
            return (a, b)

        if isinstance(b, int):
            b = (b,)

        if len(a.shape) == 1:                # 1D array is a special case
            dim = b[-1]
            while a.shape[0] > dim:          # only shrink a
                if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                    a = a[1:-1]
                else:                        # or average adjacent cells
                    a = 0.5*(a[1:] + a[:-1])
        else:
            for dim_idx in range(-(len(a.shape)),0):
                dim = b[dim_idx]
                a = a.swapaxes(0,dim_idx)        # put working dim first
                while a.shape[0] > dim:          # only shrink a
                    if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                        a = a[1:-1,:]
                    if (a.shape[0] - dim) == 1:  # or average adjacent cells
                        a = 0.5*(a[1:,:] + a[:-1,:])
                a = a.swapaxes(0,dim_idx)        # swap working dim back

        return a
        
    def readTXT(self,fname,sep=None):
        """
        Reads a txt file into an array of floats
        """
        
        fp = file(fname,'rt')
        data = np.array([map(float,line.split(sep)) for line in fp])
        fp.close()
        
        return data
    
#### For testing ####        
if __name__ == "__main__":
    BL = blend()
    BL.model_velocity()
