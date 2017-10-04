# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 11:16:08 2015

@author: dongyu
"""


import os
import subprocess
from DownloadTool import downloadUSGSriver, downloadNOAATide, downloadROMS
from FileTool import write2db, increaseT1, increaseT2, increaseT4, convertToDate, convertFormat, reviseDatFile
from suntans_driver import generateBI
from WindTool import TAMU_NCEP_wind

import pdb


def runSUNTANS(starttime, endtime, OBC_opt, IC_opt):
    """
    Specify the time for a model run and import envir var
    for example: starttime='2014-08-21-00', endtime='2014-08-23-00'
    """
    
    ## adding two days' spin up time for SUNTANS run
    starttime, endtime = increaseT1(starttime, endtime)
    
    if OBC_opt == 'ROMS'or OBC_opt == 'ROMSFILE' or IC_opt == 'ROMS':
        #download ROMS data for IC and BC
        (romst1,romst2)=increaseT4(starttime,endtime)
        downloadROMS(romst1,romst2,ROMSsrc='forecast')

    (start,end)=increaseT2(starttime,endtime)
    downloadUSGSriver(start,end)              	        #download the river inflow data
    
    staid='8772447' 					#specify the station id to download tide data 
    downloadNOAATide(start,end,staid) 	         	#download the tide
    
    #Update and write the resulting data into the database file
    dbfile = 'GalvestonObs.db'
    write2db(dbfile)
    
    """
    download the wind data from different sources and save the data in the directory DATA/ 
    """
    start1, end1 = increaseT4(starttime, endtime)
    
    print "Downloading TAMU-NCEP wind !!!\n"
    TNW = TAMU_NCEP_wind()
    TNW.writeSUNTANS('DATA/NCEP.nc')

    """
    run the shell script: buildFolder to build the rundata folder 
    and copy the necessary files to start a new run
    """
   
    subprocess.Popen('./SUNTANS/buildFolder',shell=False)

    '''
    Generate the BC and IC, 
    Time format
    start = '20141202.000000'; end = '20141206.000000'
    '''

    (t1,t2)=increaseT2(starttime,endtime)
    t1=convertFormat(t1)+'.000000'
    t2=convertFormat(t2)+'.000000'
    generateBI(t1,t2, OBC_opt, IC_opt)
    
    '''
    Revise the file suntans.dat to specify the running period and starttime
    Time format: starttime='2014-12-03', endtime='2014-12-08'
    '''
    start2, end2 = convertToDate(starttime, endtime)
    reviseDatFile(start2,end2)

    ##########################run the model############################
    dire=os.getcwd()
    os.chdir(dire+'/SUNTANS')
    #subprocess.call('make test &> suntans.out & ',shell=True)
    subprocess.call('make test ', shell=True)
    
    os.chdir(dire)
    
###########################################test only#################################################
#starttime='2014-08-21'
#endtime='2014-08-23'
#runSUNTANS(starttime, endtime)
