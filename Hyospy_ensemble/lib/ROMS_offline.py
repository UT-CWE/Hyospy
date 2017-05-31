# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:59:03 2017

@author: dongyu
"""

"""
TAMU server is unstable sometimes
This script download the data in advance of any spill incident
"""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import os
import sys

import SUNTANS
from romsio import MFncdap, roms_grid 
import othertime

import pdb

class offline_download(roms_grid):
    
    def __init__(self, starttime, endtime, **kwargs):
        
        self.__dict__.update(kwargs)
              
        
        self.t0 = datetime.strptime(starttime, '%Y-%m-%d-%H')
        self.t1 = datetime.strptime(endtime, '%Y-%m-%d-%H')  
        
        ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/NcML/oof_archive_agg', \
		'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/oof_latest_agg/roms_his_f_latest.nc']

        
        self.grdfile = '../DATA/txla_grd_v4_new.nc'    
    
        ncfile0 = [ncfiles[0]]
        ncfile1 = [ncfiles[1]]
        # Multifile object
        ftime0 = MFncdap(ncfile0,timevar='ocean_time')
        ftime1 = MFncdap(ncfile1,timevar='ocean_time')
        
        self.ind0 = othertime.findNearest(self.t0,ftime0.time)
        self.ind1 = othertime.findNearest(ftime0.time[-1], ftime0.time)
	    
        self.ind2 = othertime.findNearest(ftime0.time[-1]+timedelta(hours=1), ftime1.time)
        self.ind3 = othertime.findNearest(self.t1, ftime1.time)
	
        if (ftime1.time[self.ind2] - ftime0.time[self.ind1]).seconds/3600.>1:
	       sys.exit("ROMS forecast/hindcast data unavailable. Try again Later !!")
	
        time0 = ftime0.time[self.ind0:self.ind1+1]
        time1 = ftime1.time[self.ind2:self.ind3+1]
        self.time = np.asarray(time0.tolist()+time1.tolist())
        
        tind0, fname0 = ftime0(time0)
        tind1, fname1 = ftime1(time1)
        
        self.tind = tind0 + tind1
        self.fname = fname0 + fname1
        self.Nt = len(self.tind)
        
        ## Read the grid
        roms_grid.__init__(self, self.grdfile)        
        
        # Read the vertical coordinate variables
        self.ReadVertCoords()
                
        #pdb.set_trace()
        
    def ReadVertCoords(self):
        """
        Vertical coordinates
        """
	
        ## Vertical coordinate is from ROMS history data
        fname = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
        nc = Dataset(fname)
        

        self.Cs_r = nc.variables['Cs_r'][:]
        self.Cs_w = nc.variables['Cs_w'][:]
        self.s_rho = nc.variables['s_rho'][:]
        self.s_w = nc.variables['s_w'][:]
        self.hc = nc.variables['hc'][:]
        self.Vstretching = nc.variables['Vstretching'][:]
        self.Vtransform = nc.variables['Vtransform'][:]
        self.theta_s = nc.variables['theta_s'][:]
        self.theta_b = nc.variables['theta_b'][:]
	
        nc.close()        
        
        
    def ReadData(self, tstep):
        """
        Read variables ocean_time, zeta, temp, salt, u, v
        """                
        
        fname = self.fname[tstep]
        t0 = self.tind[tstep]
        
        print 'Reading ROMS data at time: %s...'%datetime.strftime(self.time[tstep],'%Y-%m-%d %H:%M:%S')
        
        nc = Dataset(fname)

        self.ocean_time = nc.variables['ocean_time'][t0]
        self.zeta = nc.variables['zeta'][t0,:,:]
        self.temp = nc.variables['temp'][t0,:,:,:]
        self.salt = nc.variables['salt'][t0,:,:,:]
        self.u = nc.variables['u'][t0,:,:,:]
        self.v = nc.variables['v'][t0,:,:,:]

        nc.close()


    def Writefile(self,outfile,verbose=True):
        """
        Writes subsetted grid and coordinate variables to a netcdf file
        
        Code modified from roms.py in the Octant package
        """
        self.outfile = outfile
        
        Mp, Lp = self.lon_rho.shape
        M, L = self.lon_psi.shape
        
        N = self.s_rho.shape[0] # vertical layers
        N2 = self.s_w.shape[0]
       
        xl = self.lon_rho[self.mask_rho==1.0].ptp()
        el = self.lat_rho[self.mask_rho==1.0].ptp()
        
        # Write ROMS grid to file
        nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        nc.Description = 'ROMS subsetted history file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()
        nc.type = 'ROMS HIS file'
        
        nc.createDimension('xi_rho', Lp)
        nc.createDimension('xi_u', L)
        nc.createDimension('xi_v', Lp)
        nc.createDimension('xi_psi', L)
        
        nc.createDimension('eta_rho', Mp)
        nc.createDimension('eta_u', Mp)
        nc.createDimension('eta_v', M)
        nc.createDimension('eta_psi', M)
        
        nc.createDimension('s_rho', N)
        nc.createDimension('s_w', N2)
        nc.createDimension('ocean_time', None)      
        
        nc.createVariable('xl', 'f8', ())
        nc.variables['xl'].units = 'meters'
        nc.variables['xl'] = xl
        
        nc.createVariable('el', 'f8', ())
        nc.variables['el'].units = 'meters'
        nc.variables['el'] = el
        
        nc.createVariable('spherical', 'S1', ())
        nc.variables['spherical'] = 'F'
        
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
        write_nc_var(self.angle, 'angle', ('eta_rho', 'xi_rho'))
        write_nc_var(self.h, 'h', ('eta_rho', 'xi_rho'), 'meters')
        write_nc_var(self.pm, 'pm', ('eta_rho', 'xi_rho'))
        write_nc_var(self.pn, 'pn', ('eta_rho', 'xi_rho'))
        
        write_nc_var(self.mask_rho, 'mask_rho', ('eta_rho', 'xi_rho'))
        write_nc_var(self.mask_u, 'mask_u', ('eta_u', 'xi_u'))
        write_nc_var(self.mask_v, 'mask_v', ('eta_v', 'xi_v'))
        write_nc_var(self.mask_psi, 'mask_psi', ('eta_psi', 'xi_psi'))
        
        write_nc_var(self.lon_rho, 'lon_rho', ('eta_rho', 'xi_rho'), 'meters')
        write_nc_var(self.lat_rho, 'lat_rho', ('eta_rho', 'xi_rho'), 'meters')
        write_nc_var(self.lon_u, 'lon_u', ('eta_u', 'xi_u'), 'meters')
        write_nc_var(self.lat_u, 'lat_u', ('eta_u', 'xi_u'), 'meters')
        write_nc_var(self.lon_v, 'lon_v', ('eta_v', 'xi_v'), 'meters')
        write_nc_var(self.lat_v, 'lat_v', ('eta_v', 'xi_v'), 'meters')
        write_nc_var(self.lon_psi, 'lon_psi', ('eta_psi', 'xi_psi'), 'meters')
        write_nc_var(self.lat_psi, 'lat_psi', ('eta_psi', 'xi_psi'), 'meters')
        
        # Vertical coordinate variables
        write_nc_var(self.s_rho, 's_rho', ('s_rho'))
        write_nc_var(self.Cs_r, 'Cs_r', ('s_rho'))

        write_nc_var(self.s_w, 's_w', ('s_w'))
        write_nc_var(self.Cs_w, 'Cs_w', ('s_w'))

        write_nc_var(self.hc, 'hc', ())
        write_nc_var(self.Vstretching, 'Vstretching', ())
        write_nc_var(self.Vtransform, 'Vtransform', ())
        write_nc_var(self.theta_s, 'theta_s', ())
        write_nc_var(self.theta_b, 'theta_b', ())

        #write_nc_var(self.s_rho, 's_rho', ('s_rho'))
        #write_nc_var(self.Cs_r, 'Cs_r', ('s_rho'))
        
        # Create the data variables
        create_nc_var('ocean_time',('ocean_time'),'seconds since 1970-01-01 00:00:00')
        create_nc_var('zeta',('ocean_time','eta_rho','xi_rho'),'meter')
        create_nc_var('salt',('ocean_time','s_rho','eta_rho','xi_rho'),'psu')
        create_nc_var('temp',('ocean_time','s_rho','eta_rho','xi_rho'),'degrees C')
        create_nc_var('u',('ocean_time','s_rho','eta_u','xi_u'),'meter second-1')
        create_nc_var('v',('ocean_time','s_rho','eta_v','xi_v'),'meter second-1')
        
        nc.close()    
        
    def Writedata(self, tstep):
        
        nc = Dataset(self.outfile, 'a')
        
        nc.variables['ocean_time'][tstep]=self.ocean_time
        nc.variables['zeta'][tstep,:,:]=self.zeta
        nc.variables['salt'][tstep,:,:,:]=self.salt
        nc.variables['temp'][tstep,:,:,:]=self.temp
        nc.variables['u'][tstep,:,:,:]=self.u
        nc.variables['v'][tstep,:,:,:]=self.v
	
        nc.close()
        
    def Go(self):
        """
        Downloads and append each time step to a file
        """
        if self.Nt == 0:
            sys.exit("Date and ROMS file don't match, change the input date and try again !!")
	
        for ii in range(0,self.Nt):
            self.ReadData(ii)
            self.Writedata(ii)
            
        print '##################\nDone!\n##################'
        
                



#### For testing only
if __name__ == "__main__":
    starttime='2017-05-29-00'
    endtime='2017-06-02-00'
    #offline_download(starttime, endtime)
    roms = offline_download(starttime, endtime)
    outfile = '../DATA/ROMS_offline_data.nc'
    roms.Writefile(outfile)
    roms.Go()    
    


pdb.set_trace()
