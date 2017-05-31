# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:50:14 2017

@author: dongyu
"""

"""
TAMU server is unstable when they are updating the data. 
This script continuously downloads ROMS data from their server when the data is accessible
    and dump the data in the aggregation. 
    
Data is saved in terms of month
starttime is 4 years ago, endtime is now 
Data is separated as hindcast data and forecast data
"""
import datetime
import numpy as np
from netCDF4 import Dataset
import os

import pdb

timeformat = '%Y-%m-%d-%H'
starttime = datetime.datetime.strptime('2016-01-01-00', timeformat)
end_tem = datetime.datetime.now()
endtime = datetime.datetime.strftime(end_tem, timeformat)
endtime = datetime.datetime.strptime(endtime, timeformat)

diff_hours = (endtime-starttime).total_seconds()/3600
date_list = [starttime + datetime.timedelta(hours=x) for x in range(0, int(diff_hours))]
date_list = np.asarray(date_list)

def ToMonthDate(d):
    d_by_month = datetime.datetime(d.year, d.month, 1)    
    return d_by_month
    
dates_by_month = np.array(map(ToMonthDate, date_list))
keys2 = np.unique(dates_by_month)

groups = []
for key in keys2:
    groups.append(date_list[dates_by_month==key])

for group in groups:
    start = group[0]
    end = group[-1]
    pdb.set_trace()
    filename = 'TXLA_HIS_%d%02d.nc'%(start.year, start.month)
    
    ncfiles = ['http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg']
    bbox = [-98.96,-87.80,26.00,32.24]
    ## check if the file was downloaded
    if os.path.isfile(filename):
        nc = Dataset(filename, 'r')
        
    
    
    
    

pdb.set_trace()




