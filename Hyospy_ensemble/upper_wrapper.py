# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:25:17 2017

@author: dongyu
"""
import os
import utm
import glob
import shutil
import time
import lib
import lib.SUNTANS
import hydro_wrapper
import oilspill_wrapper
from DownloadTool import downloadROMS
from FileTool import increaseT2
from WindTool import NCEP_wind, TAMU_NCEP_wind
from Probability_map import Pmap
from Blended_model import blend
import pdb


class upper_wrapper(object):
    """
    The uppermost level of wrapper
    """

    # Switches to different models
    hydro_model='BLENDED' # hydrodynamic model options: 'SUNTANS', 'ROMS', 'BLENDED'

    # ROMS data source
    ROMS_datasource='online'

    # General options
    ishindcast=True

    runGNOME=False
    runTracPy=False

    ## GNOME settings
    gnome_subset=False
    gnome_bbox=None

    probability_map=False
    google_earth=False
    mpl = 8

    number = 10 # number of ensembles
    interval = 10800 # seconds
    
    OBC_opt = 'file'  # Type 3 boundary condition option: 'constant',
    #'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'
    IC_opt = 'SUNTANS'
    
    def __init__(self,**kwargs):
        
        self.__dict__.update(kwargs)
    
    def __call__(self,starttime, endtime, period, init_latlon=[28.353786, -95.315109]):

        self.starttime = starttime
        self.endtime = endtime
        self.period = period
        self.init_latlon = init_latlon
    	
	
        if self.hydro_model == 'SUNTANS':
            self.run_suntans()

        elif self.hydro_model == 'ROMS':
            self.run_roms()

        elif self.hydro_model == 'BLENDED':
            self.run_blended()
        
        else:
            raise Exception, 'Need to set the hydro_model parameter !!'
	
    def run_suntans(self):
        """
        use SUNTANS velocity to run GNOME
        """
        
        (utm_x,utm_y)=utm.from_latlon(self.init_latlon[0], self.init_latlon[1])[0:2]
        
        for i in range(self.number):
            
            start = time.time()         
            #### run SUNTANS ####
#            hydro_wrapper.runSUNTANS(self.starttime, self.endtime, self.OBC_opt, self.IC_opt, ROMS_datasource=self.ROMS_datasource)
	
            ## Collect SUNTANS file
            basedir = os.getcwd()
            os.chdir(basedir+'/SUNTANS/rundata')
            ncfiles = []
            for ff in glob.glob("GalvCoarse_0*"):
                ncfiles.append(ff)
            os.chdir(basedir)
            
            SUNTANS_file = []
            for f in ncfiles:
                SUNTANS_file.append('%s/%s'%(basedir+'/SUNTANS/rundata', f))	
            
            ## Prepare for GNOME run
            GNOME_dir = "GNOME/%s"%str(i)   
            if os.path.exists(GNOME_dir):
                shutil.rmtree(GNOME_dir)
            os.makedirs(GNOME_dir)
                
            SUNTANS_out = '%s/txsuntans.nc'%GNOME_dir
            oilspill_wrapper.init_model(i, opt='SUNTANS')
            oilspill_wrapper.Tx_SUNTANS(SUNTANS_file, SUNTANS_out)
        	
            print "Forecast simulation, downloading TAMU-NCEP wind !!!\n"
            subset_wind = True
            TNW = TAMU_NCEP_wind(subset_wind)
            TNW.writeGNOME('%s/wind.nc'%GNOME_dir) 
            
 
            #### run GNOME ####
            if self.runGNOME:
                print 'running GNOME !!!\n'
                oilspill_wrapper.run_mul_GNOME(i, utm_x, utm_y, self.starttime, self.period, 900, opt='SUNTANS')
                oilspill_wrapper.GNOME_GM_visualization(i, opt='SUNTANS')
                oilspill_wrapper.GNOME_GE_animation(i, self.starttime, opt='SUNTANS')
                
            end = time.time()
            simulation_time = int(end-start)            
            #### pause a while for new data to come in ####
            if i != self.number-1:
                sleep_time = self.interval - simulation_time
                for j in xrange(int(sleep_time/60.)*60, 0, -60):
                    mm = j/60
                    print 'Starting new simulation in %d minutes ...'%mm
                    time.sleep(60)            
                #time.sleep(self.interval-simulation_time)
                         
        if self.probability_map:
            oilspill_wrapper.ensemble_combination(self.number, opt='SUNTANS')
            print 'creating probability map!!!\n'
            bbox=[-95.22,-94.44,28.80,29.85] # the map range
            Pmap('GNOME_combined.nc', 400, 400, self.starttime, bbox, self.mpl, self.google_earth)
	    
    

    def run_roms(self):
        """
        use ROMS velocity to run GNOME 
        A good testing date that the particles will hit SUNTANS domain is 2014-08-22~2014-08-29
        A good testing initial location is 28.353786, -95.315109	
        """
        (utm_x,utm_y)=utm.from_latlon(self.init_latlon[0], self.init_latlon[1])[0:2]        

        for i in range(self.number):
            
            start = time.time()
            #### download ROMS ####
            print "Running simulation #%s !!\n"%str(i)
#            downloadROMS(self.starttime, self.endtime, self.ROMS_datasource, ROMSsrc='forecast')
            
            ## Prepare for GNOME run
            GNOME_dir = "GNOME/%s"%str(i)   
            if os.path.exists(GNOME_dir):
                shutil.rmtree(GNOME_dir)
            os.makedirs(GNOME_dir)
            
            ROMS_file='DATA/txla_subset_HIS.nc'    
            ROMS_out = '%s/hiroms_ss_rho.nc'%GNOME_dir
            oilspill_wrapper.init_model(i, opt='ROMS')
            oilspill_wrapper.HIROMS(ROMS_file, ROMS_out) 
            
            #### wind ####
            print "Forecast simulation, downloading TAMU-NCEP wind !!!\n"
            subset_wind = False
            TNW = TAMU_NCEP_wind(subset_wind)
            TNW.writeGNOME('%s/wind.nc'%GNOME_dir)
            
            #### run GNOME ####
            if self.runGNOME:
                print 'running GNOME !!!\n'
                oilspill_wrapper.run_mul_GNOME(i, utm_x, utm_y, self.starttime, self.period, 900, opt='ROMS')
                oilspill_wrapper.GNOME_GM_visualization(i, opt='ROMS')
                oilspill_wrapper.GNOME_GE_animation(i, self.starttime, opt='ROMS')
        
            #### run TracPy ####
            if self.runTracPy:
                print 'running TracPy !!!\n'
                oilspill_wrapper.TRACPY(utm_x, utm_y, self.starttime, self.period, opt='ROMS')
            
            ## timer
            end = time.time()
            simulation_time = int(end-start)
            #### pause a while for new data to come in ####
            if i != self.number-1:
                sleep_time = self.interval - simulation_time
                for j in xrange(int(sleep_time/60.)*60, 0, -60):
                    mm = j/60
                    print 'Starting new simulation in %d minutes ...'%mm
                    time.sleep(60)
                
        
        #### probability map ####
        if self.runGNOME and self.probability_map:
            oilspill_wrapper.ensemble_combination(self.number, opt='ROMS')
            print 'creating probability map!!!\n'
            bbox=[-95.97,-94.025,27.24,29.89] # the map range
            Pmap('GNOME_combined.nc', 400, 400, self.starttime, bbox, self.mpl, self.google_earth)
            

   
    def run_blended(self):
        """
        use blended model velocity to run GNOME
        """
        (utm_x,utm_y)=utm.from_latlon(self.init_latlon[0], self.init_latlon[1])[0:2]
        
        for i in range(self.number):
            
            start = time.time()
            print "Running simulation #%s !!\n"%str(i)

            ## Step One: run SUNTANS
#            hydro_wrapper.runSUNTANS(self.starttime, self.endtime, 'ROMSFILE', 'ROMS', ROMS_datasource=self.ROMS_datasource)
	
            ## Step Two: Blend SUNTANS and ROMS
#            BL = blend(self.starttime, self.endtime)
#            BL.model_velocity()
            
            ## Prepare for GNOME run
            GNOME_dir = "GNOME/%s"%str(i)   
            if os.path.exists(GNOME_dir):
                shutil.rmtree(GNOME_dir)
            os.makedirs(GNOME_dir)

            blended_file = 'DATA/blended_uv.nc'
            blended_out = '%s/hiroms_ss_rho.nc'%GNOME_dir
            oilspill_wrapper.init_model(i, opt='blended')
            oilspill_wrapper.HIROMS(blended_file, blended_out, subset=self.gnome_subset, bbox=self.gnome_bbox)
            
            ## GNOME wind 
            print "Forecast simulation, downloading TAMU-NCEP wind !!!\n"
            subset_wind = False
            TNW = TAMU_NCEP_wind(subset_wind)
            TNW.writeGNOME('%s/wind.nc'%GNOME_dir)
            
            ## Step Three: run GNOME
            if self.runGNOME:        
                print 'running GNOME !!!\n'
                oilspill_wrapper.run_mul_GNOME(i, utm_x, utm_y, self.starttime, self.period, 900, opt='blended')
                oilspill_wrapper.GNOME_GM_visualization(i, opt='blended')
                oilspill_wrapper.GNOME_GE_animation(i, self.starttime, opt='blended')
          
            #### run TracPy ####
            if self.runTracPy:
                print 'running TracPy !!!\n'
                oilspill_wrapper.TRACPY(utm_x, utm_y, self.starttime, self.period, opt='blended')
        
            ## timer
            end = time.time()
            simulation_time = int(end-start)
            #### pause a while for new data to come in ####
            if i != self.number-1:
                sleep_time = self.interval - simulation_time
                for j in xrange(int(sleep_time/60.)*60, 0, -60):
                    mm = j/60
                    print 'Starting new simulation in %d minutes ...'%mm
                    time.sleep(60)        
        

        if self.probability_map:
            oilspill_wrapper.ensemble_combination(self.number, opt='blended')
            print 'creating probability map!!!\n'
            bbox=[-95.97,-94.025,27.24,29.89] # the map range
            Pmap('GNOME_combined.nc', 400, 400, self.starttime, bbox, self.mpl, self.google_earth)
	    
    
def timer(start, end, number, interval):
    """
    pause time
    """
    simulation_time = int(end-start)
    #### Start pausing ####
    if i!= number-1:
        sleep_time = interval - simulation_time
        for j in xrange(int(sleep_time/60.)*60, 0, -60):
            mm = j/60
            print 'Starting new simulation in %d minutes ...'%mm
            time.sleep(60) 
	




#### For testing only
if __name__ == "__main__":
    starttime='2016-03-15-00'
    endtime='2016-03-19-00'
    endtime='2016-03-16-00'
    #starttime='2017-03-15-00'
    #endtime='2017-03-19-00'
    UW = upper_wrapper()
    UW(starttime, endtime, 20, init_latlon=[28.353786, -95.315109])  #ROMS domain
    #UW(starttime, endtime, 90, init_latlon=[29.463089, -94.843460])  #SUNTANS domain
