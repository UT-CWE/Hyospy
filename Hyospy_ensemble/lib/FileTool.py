# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 14:58:33 2015

@author: dongyu
"""
from datetime import datetime, timedelta
import netcdfio
import os,sys
import shutil
from shutil import Error
from shutil import copystat
from shutil import copy2

def copy_file(src_file, dst_dir):
    """
    This function is used to move file between folders
    """
    if os.path.isfile(src_file):
        if os.path.isdir(dst_dir):
            pass
        else:
            os.makedirs(dst_dir)
        #print src_file
        srcname = src_file 
        filename = os.path.basename(src_file)
        dstname = os.path.join(dst_dir, filename)
        
        if os.path.isfile(srcname):
            copy2(srcname, dstname) 
            #print srcname,dstname,'success'
        elif os.path.isdir(dstname):
            os.remove(dstname)
            #print 'remove %s' % dstname
            copy2(srcname, dstname)
    if os.path.isfile(src_file):
        os.remove(srcname)
    if __name__ == '__main__':
        if len(sys.argv) != 3:
            print 'need srcFile and dstDir'
            sys.exit(-1)
        srcFile = sys.argv[1]
        dstDir = sys.argv[2]
        copy_file(srcFile, dstDir)

def move_file(src_file, dst_dir):
    """
    This function is used to move file between folders without deleting the old one
    """
    if os.path.isfile(src_file):
        if os.path.isdir(dst_dir):
            pass
        else:
            os.makedirs(dst_dir)
        #print src_file
        srcname = src_file 
        filename = os.path.basename(src_file)
        dstname = os.path.join(dst_dir, filename)
        
        if os.path.isfile(srcname):
            copy2(srcname, dstname) 
            #print srcname,dstname,'success'
        elif os.path.isdir(dstname):
            os.remove(dstname)
            #print 'remove %s' % dstname
            copy2(srcname, dstname)
    if __name__ == '__main__':
        if len(sys.argv) != 3:
            print 'need srcFile and dstDir'
            sys.exit(-1)
        srcFile = sys.argv[1]
        dstDir = sys.argv[2]
        copy_file(srcFile, dstDir)


"""
These two functions are used to change time to make the download data longer
than used data
Note that increaseT1 results in a time period for downloading the data
increase 2 results in a time period for generating the boundary and initial conditions
"""

def increaseT1(starttime,endtime):
    timeformat='%Y-%m-%d-%H'
    t1=datetime.strptime(starttime,timeformat)-timedelta(days=0)
    t2=datetime.strptime(endtime,timeformat)+timedelta(days=0)
    return datetime.strftime(t1, '%Y-%m-%d-%H'), datetime.strftime(t2, '%Y-%m-%d-%H')


def increaseT2(starttime,endtime):
    timeformat='%Y-%m-%d-%H'
    t1=datetime.strptime(starttime,timeformat)-timedelta(days=0)
    t2=datetime.strptime(endtime,timeformat)+timedelta(days=1)
    return (str(t1.date()),str(t2.date()))

def increaseT3(starttime,endtime):
    timeformat='%Y-%m-%d-%H'
    t1=datetime.strptime(starttime,timeformat)-timedelta(days=1)
    t2=datetime.strptime(endtime,timeformat)+timedelta(days=2)
    return (str(t1.date()),str(t2.date()))

def increaseT4(starttime,endtime):
    timeformat='%Y-%m-%d-%H'
    t1=datetime.strptime(starttime,timeformat)-timedelta(days=0)
    t2=datetime.strptime(endtime,timeformat)+timedelta(days=1)
    return datetime.strftime(t1, '%Y-%m-%d-%H'), datetime.strftime(t2, '%Y-%m-%d-%H')

def convertToDate(starttime,endtime):
    timeformat='%Y-%m-%d-%H'
    t1=datetime.strptime(starttime,timeformat) 
    t2=datetime.strptime(endtime,timeformat)
    return (str(t1.date()),str(t2.date()))


def write2db(dbfile):
    '''
    The function below aims at writing the resulting river inflow data
    and tide data into the database file
    ''' 
    #Update and write the resulting data into the database file
    create = True
    #dbfile = 'GalvestonObs.db'
    #dire1=os.path.dirname(os.path.abspath(__file__))
    dire1=os.getcwd()
    #dire1=os.path.dirname(dire2)
    folder = dire1+'/DATA/'
    #pdb.set_trace()
    ncfiles = [folder+'USGS_Rivers.nc',folder+'TCOONTide.nc']

    if create:
        print 'Creating database: %s'%dbfile
        netcdfio.createObsDB(dbfile)
    
    for nc in ncfiles:
        print 'Inserting metadata from: %s'%nc    
        netcdfio.netcdfObs2DB(nc,dbfile)

    print 'Done.'
    if os.path.isfile(folder+dbfile):
        os.remove(folder+dbfile) 
    copy_file(dbfile, folder)

"""
The following three functions are used to revise the time variable in suntans.dat
"""

def convertFormat(tt):
    return tt[0:4]+tt[5:7]+tt[8:]

def reviseDatFile(starttime,endtime):
    nt=30
    timeformat='%Y-%m-%d'
    diff=datetime.strptime(endtime,timeformat)-datetime.strptime(starttime,timeformat)
    print diff.days
    nnstep=int(diff.days)*3600*24/nt
    print nnstep
    string1='20160'
    string2='20141203.010000'
    dire=os.getcwd()
    f=open(dire+'/SUNTANS/rundata/suntans1.dat','r')
    content=f.readlines()
    open(dire+'/SUNTANS/rundata/suntans.dat','w').write('')
    h=open(dire+'/SUNTANS/rundata/suntans.dat','w')
    newtime=convertFormat(starttime) +'.010000'   ##SUNTANS starttime
    for line in content:
        h.write(line.replace(string1,str(nnstep)).replace(string2,newtime))
    h.close
'''
    for line in content:
    #print line
        if 'starttime' in line:
            line= 'starttime' + ' ' + convertFormat(starttime) +'.010000'+'    '
        content1.append(line)
    content1[37]='nsteps'+'  '+ str(nnstep) + '           '    
    for line in content1:
        open(dire+'/CoarseTri/rundata/suntans.dat','a').write(line)
'''
def reviseDatFile2(starttime,endtime):
    nt=1
    timeformat='%Y-%m-%d'
    diff=datetime.strptime(endtime,timeformat)-datetime.strptime(starttime,timeformat)
    print diff.days
    nnstep=int(diff.days)*3600*24/nt
    print nnstep
    string1='172800'
    string2='20141102.000000'
    dire=os.getcwd()
    f=open(dire+'/FineTri/rundata/suntans1.dat','r')
    content=f.readlines()
    open(dire+'/FineTri/rundata/suntans.dat','w').write('')
    h=open(dire+'/FineTri/rundata/suntans.dat','w')
    newtime=convertFormat(starttime) +'.000000'
    for line in content:
        h.write(line.replace(string1,str(nnstep)).replace(string2,newtime))
    h.close


