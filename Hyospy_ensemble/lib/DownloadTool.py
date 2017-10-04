# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 16:05:51 2015

@author: dongyu
"""
import urllib2
import string
import datetime
import numpy as np
from lxml import etree
from xml.etree import ElementTree
import os
import netcdfio
import urllib
import shutil
import tarfile
from contextlib import closing
import math
from netCDF4 import Dataset
import time as tim
from othertime import MinutesSince
from romsio import roms_subset
from FileTool import copy_file
import pdb




def downloadUSGSriver(starttime,endtime):
    from getUSGSnwis import getUSGSnwis
    stationids = ['08066500',\
                '08078000',\
                '08067500',\
                '08073600',\
                '08042558',\
                '08031000',\
                '08067525']
      
    ncfile = 'USGS_Rivers.nc'
    getUSGSnwis(stationids,starttime,endtime,ncfile)    
    files=['USGS_Rivers.dbf','USGS_Rivers.nc','USGS_Rivers.shp','USGS_Rivers.shx']
    dst_dir=os.getcwd()+'/DATA/'
    for src_file in files:
        copy_file(src_file, dst_dir)


def downloadNOAATide(starttime,endtime,staid):
    data=downloadTide(starttime,endtime,staid)
    dire=os.path.dirname(__file__)
    dire=os.getcwd()
    ncfile = dire+'/DATA/TCOONTide.nc'
    globalatts = {'title':'TCOON oceanographic observation data'}
    netcdfio.writePointData2Netcdf(ncfile,data,globalatts)
    print 'end writing tide data...'
    

def downloadTide(starttime,endtime,staid):
    """ Main function to construct data in a uniform format"""
    lon=-95.3025
    lat=28.9433
    ID=8772447
    nn='Galveston Freeport'
    vv='waterlevel'
    # Build up the output data as a list of dictionaries
    meta=[]
    coords = [{'Name':'longitude','Value':lon,'units':'degrees East'},\
        {'Name':'latitude','Value':lat,'units':'degrees North'},\
        {'Name':'time','Value':[],'units':'minutes since 1970-01-01 00:00:00'}]
    attribs = {'StationID':str(ID),'StationName':nn,'Data':[],'coordinates':'time, longitude, latitude','coords':coords} 
    meta.append(attribs)

    data=[]
    tmp = {vv:meta[0]}
    ctr=0

    output,t,atts = getNOAATide(starttime,endtime,staid)
    
    if ctr==0:                
        tmp[vv].update(atts)
    # Append the data to the list array
    tmp[vv]['Data'] += output
    # Append the time data
    ctr=-1
    for cc in tmp[vv]['coords']:
        ctr+=1
        if cc['Name']=='time':
            tmp[vv]['coords'][ctr]['Value'] += t 
    if np.size(tmp[vv]['Data']) > 0:
            data.append(tmp)
            
    return data

def getTCOON(starttime,endtime,staid,choice=True):
    startt=starttime[5:7]+'.'+starttime[8:10]+'.'+starttime[0:4]
    endt=endtime[5:7]+'.'+endtime[8:10]+'.'+endtime[0:4]
    url = 'http://lighthouse.tamucc.edu/pd?stnlist='+staid+'&serlist=pwl%2Charmwl&when='+startt+'-'+endt+'&whentz=UTC0&-action=c&unit=metric&elev=msl'
    try:
        print 'Opening: %s'%url
        f=urllib2.urlopen(url)
    except:
        raise Exception, 'cannot open url:\n%s'%url
    data1=[]
    data0=[] #Note data0 is different from data that data0 includes the possible string 'NA'
    data=[] #To store the tidal elevation and time
    attribs = {'long_name':'Water surface elevation','units':'m'}
    for s in f:
        if '#' not in s:
            line0=s.split()
            data0.append(line0)
            if 'NA' not in s:
                line=s.split()
                data1.append(line)
    if choice:
        # calculate the mean difference between pwl and harmwl with N = 30 (3 hours)        
        sum=0
        #pdb.set_trace()
        for n in range(len(data1),len(data1)-30,-1):
            if data1[n-1][1]== 'RM': data1[n-1][1] = data1[n-1][2]
            aa=string.atof(data1[n-1][1]);bb=string.atof(data1[n-1][2])
            diff=aa-bb
            sum=sum+diff
        x=sum/len(data1)
        # find the index of the first 'NA'
        column2=[]
        for line in data0:
            column2.append(line[1])
        #pdb.set_trace()
        if 'NA' in column2:
            y=column2.index('NA')
        else:
            y=len(data0)
        # add x with harmwl to creat the forecast part of pwl, set the return interval as 72 hours (72*60/6=720)
        for ii in range(len(data0)):
            if data0[ii][1]=='RM':
                data0[ii][1]=data0[ii][2]
        for i in range(y,len(data0)):
            data0[i][1]=string.atof(data0[i][2])+x*(1-i/720)
        for k in range(y):
            data0[k][1]=string.atof(data0[k][1])
        for m in range(len(data0)):
            del data0[m][2]
        for line in data0:
            data.append(line)
    else:
        for line in data0:
            if line[1]=='NA':
                line[1]=0
            del line[2]
            data.append(line)
    tt=[] #time
    for kk in range(len(data)):
        outTime=parseTime(data[kk][0])
        tt.append(outTime)
    elev=[] #water elevation
    for ii in range(len(data)):
        tem=data[ii][1]
        elev.append(tem)
    
    print 'finish downloading data from TCOON...'
    pdb.set_trace()
    return elev, tt, attribs

def getNOAATide(starttime, endtime, staid):

    attribs = {'long_name':'Water surface elevation','units':'m'}    
    timeformat = "%Y-%m-%d"
    start = datetime.datetime.strptime(starttime, timeformat)
    end = datetime.datetime.strptime(endtime, timeformat) + datetime.timedelta(days=1)
    
    start = datetime.datetime.strftime(start, "%Y%m%d")
    end = datetime.datetime.strftime(end, "%Y%m%d")
    
    url_water_level = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=%s 00:00&end_date=%s 00:00&station=%s&product=water_level&datum=msl&units=metric&time_zone=gmt&application=web_services&format=xml'%(start, end, staid)
    url_predictions = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date=%s 00:00&end_date=%s 00:00&station=%s&product=predictions&datum=msl&units=metric&time_zone=gmt&application=web_services&format=xml'%(start, end, staid)
    
    
    fp = urllib.urlopen(url_water_level)
    tree = etree.parse(fp)
    root = tree.getroot()
    fp.close()
        
    print root[1][0].attrib['t']
    
    time = []
    water_level = []
    for i in range(len(root[1])):
        time.append(datetime.datetime.strptime(root[1][i].attrib['t'], "%Y-%m-%d %H:%M"))
        water_level.append(root[1][i].attrib['v'])
        
    water_level = np.asarray(water_level)
    water_level[water_level==''] = '0'
    ele = []
    for i in range(len(water_level)):
        ele.append(string.atof(water_level[i]))
        
    time = np.asarray(time)
    ele = np.asarray(ele)
    #### predictions ####
    fp2 = urllib.urlopen(url_predictions)
    tree2 = etree.parse(fp2)
    root2 = tree2.getroot()
    fp2.close()
        
    print root2[0].attrib['t']
    
    time2 = []
    predictions = []
    for i in range(len(root2)):
        time2.append(datetime.datetime.strptime(root2[i].attrib['t'], "%Y-%m-%d %H:%M"))
        predictions.append(root2[i].attrib['v'])
        
    predictions = np.asarray(predictions)
    predictions[predictions==''] = '0'
    ele2 = []
    for i in range(len(predictions)):
        ele2.append(string.atof(predictions[i])) 

    time2 = np.asarray(time2)
    ele2 = np.asarray(ele2)    
    
    nn = len(ele)
    ele2[:nn] = ele
    
    tt=[] #time
    for kk in range(len(time2)):
        outTime=parseTime2(time2[kk])
        tt.append(outTime)
    #pdb.set_trace()
    return ele2.tolist(), tt, attribs
    
 
def parseTime(inTime):
    basetime=datetime.datetime(2014,1,1,00,00)
    yy=(int(inTime[0:4])-2014)*365
    dd=inTime[4:7]
    hh=inTime[8:10]
    mm=inTime[10:]
    delta=datetime.timedelta(days = yy+ int(dd)-1,hours=int(hh),minutes=int(mm))
    outTime=basetime+delta
    return MinutesSince(outTime)[0]  

def parseTime2(inTime):

    return MinutesSince(inTime)[0]  


def downloadROMS(starttime,endtime, ROMS_datasource, ROMSsrc='forecast'):
    """
    This function is used to download the ROMS data from http://barataria.tamu.edu:8080/thredds/catalog.html
    The data is saved as DATA/txla_subset_HIS.nc, both for IC, BC of SUNTANS and run drifter model

        starttime = '2014-08-21-00'
        endtime   = '2014-08-22-00'
        ROMSsrc = 'forecast' or 'hindcast' 

    """
    timelims = (starttime.replace("-", "")+'0000', endtime.replace("-", "")+'0000')
    
    if ROMSsrc=='forecast': 
	"""
	6 days forecast
	"""
	
	if ROMS_datasource == 'online':
	    print "sourcing ROMS data using OpenDap!!\n"
   
            grdfile = 'http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc'
            
	    ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/NcML/oof_archive_agg', \
		    'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/oof_latest_agg/roms_his_f_latest.nc']

	elif ROMS_datasource == 'offline':
	    print "sourcing ROMS data using HTTPServer!!\n"
	    if os.path.isfile('DATA/ROMS_offline_data.nc'):
       		ncfiles = ['DATA/ROMS_offline_data.nc']
	        grdfile = 'DATA/txla_grd_v4_new.nc'
	    else:
	        raise IOError('Need to download data from TAMU server first !!')


    elif ROMSsrc=='hindcast':
        """
        Two datasets
        """
	
	timeformat='%Y-%m-%d-%H'
	if datetime.datetime.strptime(starttime, timeformat) > datetime.datetime(2017,1,1) and \
		datetime.datetime.strptime(endtime, timeformat) < datetime.datetime.now():
	    ## Lastest 3 Months History
	    ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/NcML/oof_archive_agg']

	elif datetime.datetime.strptime(endtime, timeformat) < datetime.datetime(2017,1,1):
	    ## 1993.1.1~2017.1.1
	    ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg']

	else:
	    raise IOError('Switch to forecast !!')

	grdfile = 'http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc'

    else:
        raise "There is no such option in choosing ROMS data source"
    
    bbox = [-98.96,-87.80,26.00,32.24]
 
    
    roms = roms_subset(ncfiles,bbox,timelims,gridfile=grdfile)
    outfile = os.getcwd()+'/DATA/txla_subset_HIS.nc'
    roms.Writefile(outfile)
    roms.Go()
