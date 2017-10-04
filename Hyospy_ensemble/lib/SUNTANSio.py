# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 12:05:57 2016

@author: dongyu
"""

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date
from datetime import datetime, timedelta
import othertime
import os
import pdb

class SUNTANS_download(object):
    """
    General class for downloading SUNTANS data in the format of OPENDAP
    """
    gridfile = None    
    
    def __init__(self, ncfiles, timelims, **kwargs):
        self.__dict__.update(kwargs)
        
        if self.gridfile == None:
            self.gridfile = ncfiles[0]
            
        self.ncfiles = ncfiles
        
        
        # Step 1) Find the time steps
        self.t0 = datetime.strptime(timelims[0],'%Y%m%d%H%M%S')
        self.t1 = datetime.strptime(timelims[1],'%Y%m%d%H%M%S')
        
        # Multifile object
        ftime = MFncdap(ncfiles,timevar='time')
        
        ind0 = othertime.findNearest(self.t0,ftime.time)
        ind1 = othertime.findNearest(self.t1,ftime.time)
        
        self.timei = ftime.time[ind0:ind1]        
        self.tind,self.fname = ftime(self.timei) # list of time indices and corresponding files
        
        self.Nt = len(self.tind)
        # Step 2) Read the grid variables
        self.ReadGrid(self.gridfile)  
        
    
    def ReadGrid(self, grdfile):
        """
        funtion that reads the grid variables
        """
        nc =  Dataset(grdfile,'r')
        
        self.xv = nc.variables['xv'][:]
        self.yv = nc.variables['yv'][:]
        self.xp = nc.variables['xp'][:]
        self.yp = nc.variables['yp'][:]
        self.xe = nc.variables['xe'][:]
        self.ye = nc.variables['ye'][:]
        self.dz = nc.variables['dz'][:] 
        self.dv = nc.variables['dv'][:]
        self.Ac = nc.variables['Ac'][:]
        self.Nk = nc.variables['Nk'][:]
        self.face = nc.variables['face'][:]
        self.mark = nc.variables['mark'][:]
	self.cells = nc.variables['cells'][:]
        
        self.Nc = len(self.xv)
        self.Np = len(self.xp)
        self.Ne = len(self.xe)
        self.Nk = len(self.dz)
        self.numsides = self.face.shape[1]
        
    def ReadData(self, tstep):
        """
        funtion that read data at every time step
        """
        fname = self.fname[tstep]
        t0 = self.tind[tstep]
        
        print 'Reading SUNTANS data at time: %s...'%datetime.strftime(self.timei[tstep],'%Y-%m-%d %H:%M:%S')  
        nc = Dataset(fname)
        
        self.time = nc.variables['time'][t0]
        
        self.temp = nc.variables['temp'][t0,:,:]
        self.salt = nc.variables['salt'][t0,:,:]
        self.uc = nc.variables['uc'][t0,:,:]
        self.vc = nc.variables['vc'][t0,:,:]
        self.nu_v = nc.variables['nu_v'][t0,:,:]
        self.rho = nc.variables['rho'][t0,:,:]
        self.tau_x = nc.variables['tau_x'][t0,:]
        self.tau_y = nc.variables['tau_y'][t0,:]
        self.eta = nc.variables['eta'][t0,:]
        


    def Writefile(self, outfile, verbose=True):
        """
        Write the downloaded SUNTANS file into a new netcdf file
        """
        
        self.outfile = outfile
        
        # Write SUNTANS grid to file
        nc = Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
        nc.Description = 'SUNTANS subsetted history file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()
        nc.type = 'SUNTANS HIS file'
        #pdb.set_trace()
        nc.createDimension('Nc', self.Nc)
        nc.createDimension('Np', self.Np)
        nc.createDimension('Ne', self.Ne)
        nc.createDimension('Nk', self.Nk)
        nc.createDimension('numsides', self.numsides)
    
        nc.createDimension('time', None)
        
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
        write_nc_var(self.xv, 'xv', ('Nc'))
        write_nc_var(self.yv, 'yv', ('Nc'))
        write_nc_var(self.xp, 'xp', ('Np'))
        write_nc_var(self.yp, 'yp', ('Np'))
        write_nc_var(self.xe, 'xe', ('Ne'))
        write_nc_var(self.ye, 'ye', ('Ne'))
        write_nc_var(self.dz, 'dz', ('Nk'))
        write_nc_var(self.dv, 'dv', ('Nc'))
        write_nc_var(self.Ac, 'Ac', ('Nc'))
        write_nc_var(self.Nk, 'Nk', ('Nc'))
        write_nc_var(self.face, 'face', ('Nc','numsides'))
        write_nc_var(self.mark, 'mark', ('Ne'))
        write_nc_var(self.cells, 'cells', ('Nc','numsides'))
        
        
        # Create the data variables
        create_nc_var('time',('time'),'seconds since 1990-01-01 00:00:00')
        create_nc_var('salt',('time','Nk','Nc'),'psu')
        create_nc_var('temp',('time','Nk','Nc'),'degrees C')
        create_nc_var('uc',('time','Nk','Nc'),'meter second-1')
        create_nc_var('vc',('time','Nk','Nc'),'meter second-1')
        create_nc_var('nu_v',('time','Nk','Nc'),'m2 s-1')
        create_nc_var('rho',('time','Nk','Nc'),'kg m-3')
        create_nc_var('tau_x',('time','Nc'),'N m-2')
        create_nc_var('tau_y',('time','Nc'),'N m-2')
        create_nc_var('eta',('time','Nc'),'m')
        
        nc.close()
        
    def Writedata(self, tstep):
        """
        funtion that write the data into the variables
        """
        
        nc = Dataset(self.outfile, 'a')
        
        nc.variables['time'][tstep] = self.time
        nc.variables['salt'][tstep] = self.salt
        nc.variables['temp'][tstep] = self.temp
        nc.variables['uc'][tstep] = self.uc
        nc.variables['vc'][tstep] = self.vc
        nc.variables['nu_v'][tstep] = self.nu_v
        nc.variables['rho'][tstep] = self.rho
        nc.variables['tau_x'][tstep] = self.tau_x
        nc.variables['tau_y'][tstep] = self.tau_y
        nc.variables['eta'][tstep] = self.eta
        
        nc.close()
        
    def Go(self):
        """
        Downloads and append each time step to a file
        """
        for ii in range(self.Nt):
            self.ReadData(ii)
            self.Writedata(ii)
            
        print '##################\nDone!\n##################'
        
        
        
        
        
class MFncdap(object):
    """
    Multi-file class for opendap netcdf files
    
    MFDataset module is not compatible with opendap data 
    """
    
    timevar = 'time'
    
    def __init__(self,ncfilelist,**kwargs):
        
        self.__dict__.update(kwargs)
        
        self.timelookup = {}
        self.time = np.zeros((0,))
        for f in ncfilelist:
            print f
            nc = Dataset(f)
            t = nc.variables[self.timevar]
            time = num2date(t[:],t.units)
            nc.close()
            
            self.timelookup.update({f:time})
            self.time = np.hstack((self.time,np.asarray(time)))
            
        self.time = np.asarray(self.time)
    
            
    def __call__(self,time):
        """
        Return the filenames and time index of the closest time
        """
        
        fname = []
        tind =[]
        for t in time:
            flag=1
            for f in self.timelookup.keys():

                if t >= self.timelookup[f][0] and t<=self.timelookup[f][-1]:
#                    print 'Found tstep %s'%datetime.strptime(t,'%Y-%m-%d %H:%M:%S')
                    tind.append(othertime.findNearest(t,self.timelookup[f][:]))
                    fname.append(f)
                    flag=0

#            if flag:
#                print 'Warning - could not find matching file for time:%s'%datetime.strptime(t,'%Y-%m-%d %H:%M:%S')
#                tind.append(-1)
#                fname.append(-1)
        
        return tind, fname
        
########################### Direct calling class #############################
if __name__ == "__main__":
    """
    download SUNTANS output from the server
    """
    def _int2str(num):
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
            
    starttime='2009-04-01'
    endtime='2009-05-01'
    timelims = (starttime.replace("-", "")+'000000', endtime.replace("-", "")+'000000')
    
    ncfiles = []
    
    basedir='http://barataria.tamu.edu:8080/thredds/dodsC/2009'
    #basedir = 'http://barataria.tamu.edu:8080/thredds/dodsC/FineTri_2009'
    for i in range(15,65):
    #for i in range(181, 246):
        filestr = '%s/GalvCoarse_2009_AVG_0%s.nc'%(basedir, _int2str(i))
        #filestr = '%s/GalvTri_2009_AVG_0%s.nc'%(basedir, _int2str(i))
        ncfiles.append(filestr)
        
    sun = SUNTANS_download(ncfiles,timelims)
    #outfile = os.getcwd()+'/SUNTANS_file/GalvCoarse.nc'
    #outfile = os.getcwd()+'/SUNTANS_file/GalvTri.nc'
    outfile = os.getcwd()+'/GalvCoarse200901.nc'
    sun.Writefile(outfile)
    sun.Go()
