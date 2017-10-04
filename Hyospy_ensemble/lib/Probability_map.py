# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 16:05:59 2016

@author: dongyu
https://github.com/nasa/PyAMPR/blob/master/pyampr/google_earth_tools.py
https://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/
https://github.com/hetland/octant/blob/master/octant/sandbox/googleearth.py
"""

from netCDF4 import Dataset, MFDataset
import matplotlib.pyplot as plt
import numpy as np
import shutil
import subprocess
import os
import time
from datetime import datetime
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)

import FileTool

import pdb


class Pmap(object):
    """
    A general class for visualizing the concentration of oil spill in Google Earth
    """
    def __init__(self,filename, row, col, starttime, bbox=[-95.25,-94.45,28.9,29.9], mpl=8, google_earth=True):
        """
        filename is the name of GNOME output
        row and col are the rows and columns of the grid to calculate the particle concentration
        bbox is the domain to visualize the data
        mpl is the multiplier and the interval to generate figures, default setting is 8
        for example, there are 200 files avaiable, but you need to generate 25 figures to save GPU storage and time
        """
        # change directory to GNOME
        basedir = os.getcwd()
    	os.chdir(basedir+'/GNOME/')
        
        dire=os.getcwd()
        nc=Dataset(filename,'r')      
        particle_num=int(nc.variables['particle_count'][:][0])
        data_length=len(nc.variables['longitude'][:])
        
        self.runC(filename, data_length, particle_num, bbox, row, col)
	print 'Calculating probability, this takes about 1 min !!!\n'
	for i in xrange(80, 0, -10):
	    print '%d seconds left ...'%i
	    time.sleep(10)
	    
        os.chdir(dire)
        ############read count matrix (particle trajectories)############
        west=bbox[0]; east=bbox[1]
        south=bbox[2]; north=bbox[3]
        
        x = np.linspace(west, east, col)
        y = np.linspace(south, north, row)
        X, Y = np.meshgrid(x, y)
        ############count time length AKA files number in /data_C
        DIR = dire+'/data_C'
        ntime=len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])-8
        
        
        if google_earth:
	    print "Creating probability map on Google Earth !!!\n"
	    ##########convert coordinate to google earth projection########
	    if os.path.exists(dire+'/figure'):
	        shutil.rmtree(dire+'/figure')
            os.makedirs(dire+'/figure')

            pixels = 512 * 10
            for j in range(ntime/mpl):
                i = j*mpl
                filename=dire+'/data_C/'+'out_count'+str(i)+'.nc'
                a=Dataset(filename,'r')
                count=a.variables['data'][:]
                count=count.astype('float')
                count[count==0]=np.nan      
		## calculate probability
		count = count / count[np.isnan(count)==False].max()      
    
                fig, ax = self.gearth_fig(llcrnrlon=west,
                                     llcrnrlat=south,
                                     urcrnrlon=east,
                                     urcrnrlat=north,
                                     pixels=pixels)
                     
                cs = ax.contourf(X,Y,count,200)
                ax.set_axis_off()
                plt.savefig(dire+'/figure/'+str(i)+'.png', transparent=True)
            	plt.close()

	    ## colorbar
	    fig1 = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
            ax1 = fig1.add_axes([0.0, 0.05, 0.2, 0.9])
            cb = fig1.colorbar(cs, cax=ax1)
            #cb.set_label('Probability', rotation=-90, color='k') #, labelpad=20)
            cb.set_label(' ', rotation=-90, color='k')
            fig1.savefig('legend.png', transparent=False, format='png')

            bbox=[west, south, east, north]
            #################Generate the KML file from figures#############
            self.make_kml(ntime,bbox,starttime,mpl)

            ############# Deleting figures and count_out files ############# 
            fig_folder=dire+'/figure'
            for the_file in os.listdir(fig_folder):
                file_path = os.path.join(fig_folder, the_file)
                if os.path.isfile(file_path):
                    os.unlink(file_path)
	
	else:
	    print "Creating probability map on Python Basemap!!!\n"
 	    from mpl_toolkits.basemap import Basemap
            count0 = np.zeros([row, col])
	    for j in range(ntime/mpl):
                i = j*mpl
	        filename = '%s/data_C/out_count%s.nc'%(dire, str(i))
	    	a = Dataset(filename,'r')
            	count=a.variables['data'][:]
            	count=count.astype('float')
            	#count[count==0]=np.nan
                ## calculate probability
		count = count / count.max()  
                count0 += count  
            
	    count0[count0==0]=np.nan
	    count0 = count0 / count0[np.isnan(count0)==False].max()
	    pdb.set_trace()
	    ## plotting
            fig = plt.figure(figsize=(10,8))
	    basemap = Basemap(projection='merc',llcrnrlat=south,urcrnrlat=north,\
          	llcrnrlon=west,urcrnrlon=east, resolution='h')
          
	    basemap.drawcoastlines()
   	    basemap.fillcontinents(color='coral',lake_color='aqua')
	    basemap.drawcountries()
	    basemap.drawstates()  

	    llons, llats=basemap(*(X,Y)) 
	    con=basemap.contourf(llons,llats,count0,200,zorder=4, alpha=0.6)
	    cbar = plt.colorbar(con, orientation='vertical')
	    cbar.set_label("probability")
	    pdb.set_trace()

        
        for i in range(ntime):
            filename=dire+'/data_C/'+'out_count'+str(i)+'.nc'
            self.deleteFile(filename)  

	## change directory back 
	os.chdir(basedir)

        
    def deleteFile(self,srcname):
        if os.path.isfile(srcname):
            os.remove(srcname)


    def runC(self, filename,  data_length, particle_num, bbox, row_num, column_num):
        """
        The function is used to update the C++ program with parameters 
        and execuate the program to generate netcdf file for particle concentration
        """
	

        dire=os.getcwd()
	if os.path.exists(dire+'/data_C'):
	    shutil.rmtree(dire+'/data_C')
        os.makedirs(dire+'/data_C')

	scripts = ['../DATA/data_C/Makefile', '../DATA/data_C/Makefile.in', \
			'../DATA/data_C/old_readnetcdf.cpp', filename]
	for ss in scripts:
            FileTool.move_file(ss, dire+'/data_C')

        os.chdir(dire+'/data_C')
        f=open(dire+'/data_C/old_readnetcdf.cpp','r')
        content=f.readlines()
        open(dire+'/data_C/readnetcdf.cpp','w').write('')
        h=open(dire+'/data_C/readnetcdf.cpp','w')
        for line in content:
            h.write(line.replace('data_length',str(data_length)).replace('particle_num',str(particle_num)). replace('west_point',str(bbox[0])).replace('east_point',str(bbox[1])).replace('south_point',str(bbox[2])).replace('north_point',str(bbox[3])).replace('GNOME_file',filename).replace('row_num',str(row_num)).replace('column_num',str(column_num)))
        h.close()
        subprocess.Popen('make',shell=False)
        time.sleep(3)
        subprocess.Popen('./readnetcdf',shell=False)
        
    def make_kml(self, ntime, bbox, starttime, mpl, colorbar=True):
        """TODO: LatLon bbox, list of figs, optional colorbar figure,
        and several simplekml kw..."""
        """bbox is the cooridnate"""
        llcrnrlon=bbox[0]    #west 
        llcrnrlat=bbox[1]    #south 
        urcrnrlon=bbox[2]    #east 
        urcrnrlat=bbox[3]    #north
        dire=os.getcwd()+'/figure/'
        obj=[str(i) for i in xrange(ntime)]
        kml = Kml()
        altitude = 161000               #This is 100 miles          #set the view altitude: 2e7
        roll = 0
        tilt = 0
        altitudemode = AltitudeMode.relativetoground
        camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                        longitude=np.mean([urcrnrlon, llcrnrlon]),
                        altitude=altitude, roll=roll, tilt=tilt,
                        altitudemode=altitudemode)

        kml.document.camera = camera
        
	starttime = datetime.strptime(starttime, '%Y-%m-%d-%H')
	yr = datetime.strftime(starttime, '%Y')
        month = datetime.strftime(starttime, '%m')
        day = datetime.strftime(starttime, '%d')
	hr = datetime.strftime(starttime, '%H')
        #yr=tbox[0]         #2014
        #month=tbox[1]      #8
        #day=tbox[2]        #21
        for j in range(ntime/mpl):
            i=j*mpl
            fig=dire+str(i)+'.png'            
            #draworder += 1
            obj[i] = kml.newgroundoverlay(name='GroundOverlay')
            #obj[i].draworder = draworder
            obj[i].visibility = 1
            obj[i].name = 'overlay'
            obj[i].color = '9effffff'
            obj[i].latlonbox.rotation = 0
            obj[i].gxaltitudemode = 'clampToSeaFloor'
            #obj[i].timespan.begin="2011-02-20T"+str("%02d" % i)+":00:00Z"
            #obj[i].timespan.end="2011-02-20T"+str("%02d" % (i+1))+":00:00Z"
            obj[i].timestamp.when=yr+'-'+month+'-'+str('%02d' % (int(day)+i/96))+'T'\
                                   + str("%02d" % (int(hr)+ (i/4)%24))+':' +str("%02d" %((i*15)%60))+':00Z'
                               
            obj[i].icon.href = fig
            obj[i].latlonbox.east = llcrnrlon
            obj[i].latlonbox.south = llcrnrlat
            obj[i].latlonbox.north = urcrnrlat
            obj[i].latlonbox.west = urcrnrlon
            #pdb.set_trace()

	    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
                screen = kml.newscreenoverlay(name='ScreenOverlay')
                screen.icon.href = 'legend.png'
                screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
                screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
                screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                    xunits=Units.fraction,
                                    yunits=Units.fraction)
                screen.size.x = 0
                screen.size.y = 0
                screen.size.xunits = Units.fraction
                screen.size.yunits = Units.fraction
                screen.visibility = 1
    
        kmzfile = 'GNOME_PM.kmz'
        #kmzfile = kw.pop('kmzfile', 'GNOME_GE.kmz')
        kml.savekmz(kmzfile)
        print "#### Successfully generate the Google Earth file %s"%kmzfile
       
    
    def gearth_fig(self, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
        """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
        aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
        xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
        ysize = np.ptp([urcrnrlat, llcrnrlat])
        aspect = ysize / xsize
    
        if aspect > 1.0:
            figsize = (10.0 / aspect, 10.0)
        else:
            figsize = (10.0, 10.0 * aspect)

        if False:
            plt.ioff()  # Make `True` to prevent the KML components from poping-up.
        fig = plt.figure(figsize=figsize,
                         facecolor=None,
                         frameon=False,
                         dpi=pixels//10)
        # KML friendly image.  If using basemap try: `fix_aspect=False`.
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(llcrnrlon, urcrnrlon)
        ax.set_ylim(llcrnrlat, urcrnrlat)
        return fig, ax
        
if __name__ == "__main__":
    bbox=[-95.78,-94.45,28.3,29.9]
    google_plot('GNOME_output.nc', 300, 300, bbox, mpl=8)
