# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 12:21:56 2015

"""

from netCDF4 import Dataset, num2date, date2num, date2index, MFDataset
import glob, os
import math
import datetime as dt
import numpy as np

def get_var_map(filename, var_list=['lon','lat','time','u','v','air_u','air_v']):
    
    var_map = dict()
    ncvars = Dataset(filename).variables
    long_names = { 'lon':   'longitude',
                   'lat':   'latitude',
                   'time':  'time',
                   'u':     'eastward_sea_water_velocity',
                   'v':     'northward_sea_water_velocity',
                   'air_u': 'eastward_wind',
                   'air_v': 'northward_wind',
                 }
                 
    for var in var_list:
        matches = []
        for varname in ncvars:
            try:
                if ncvars[varname].standard_name == long_names[var]:
                    matches.append(varname)
            except AttributeError:
                pass
        if not matches:
            var_map[var] = None
        elif len(matches) == 1:
            var_map[var] = matches[0]
        else:
            var_map[var] = matches
            
    return var_map

def fix_time_units(units):
    '''
    GNOME doesn't support units of this form: 'hours since 2014-12-12T18:00:00Z'
    Just replace T with space
    '''
    new_units = ' '.join(units.split('T'))
    return new_units

def show_ncfile_tbounds(filename,tvar='time'):
    
    t = Dataset(filename).variables[tvar]
    print 'Start date: ', num2date(t[0],t.units)
    try:
        print 'End date: ', num2date(t[-1],t.units)
    except IndexError:
        print num2date(t[0],t.units)

def show_tbounds(t):
    
    print 'Start date: ', num2date(t[0],t.units)
    try:
        print 'End date: ', num2date(t[-1],t.units)
    except IndexError:
        print num2date(t[0],t.units)
        
def get_tindex(t,start_date,end_date,stride=None):
        
    tindex = []
    tindex.append(date2index(start_date,t,select='before'))
    tindex.append(date2index(end_date,t,select='after') + 1)
    if stride is None:
        tindex.append(1)
    else:
        tindex.append(stride)
    return tindex

def adjust_time(t,t_units):
    #GNOME can't handle pre-1970 date units
    
    dtime = num2date(t,t_units)
    new_units = 'days since 1980-1-1 00:00:00'
    new_time = date2num(dtime,units=new_units)
    
    return new_time,new_units

def round_time(datetime_in=None, roundto=60):
   """Round a datetime object to any time laps in seconds
   dt : datetime.datetime object
   roundTo : Closest number of seconds to round to, default 1 minute.
   Author: Thierry Husson 2012 - Use it as you want but don't blame me.
   """
   seconds = (datetime_in - datetime_in.min).seconds
   # // is a floor division, not a comment on following line:
   rounding = (seconds+roundto/2) // roundto * roundto
   return datetime_in + dt.timedelta(0,rounding-seconds,-datetime_in.microsecond)
   
def make_filelist_for_GNOME(file_dir,file_match='*.*',outfilename='filelist.txt'):
    #used to load multiple files as one mover in GNOME
    flist = glob.glob(os.path.join(file_dir,file_match))
    f = open(os.path.join(file_dir,outfilename),'w')
    f.write('NetCDF Files\n')
    f.write('\n'.join(['[FILE] ' + os.path.split(file)[-1] for file in flist]))
    f.close()


def utmToLatLng(zone, easting, northing, northernHemisphere=True):
    # Convert UTM coordinates to lat/lon
    #TODO: this method probably belongs in a more generic module

    if not northernHemisphere:
        northing = 10000000 - northing

    a = 6378137
    e = 0.081819191
    e1sq = 0.006739497
    k0 = 0.9996

    arc = northing / k0
    mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 * math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))

    ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / (1 + math.pow((1 - e * e), (1 / 2.0)))

    ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0

    cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
    cc = 151 * math.pow(ei, 3) / 96
    cd = 1097 * math.pow(ei, 4) / 512
    phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + cc * math.sin(6 * mu) + cd * math.sin(8 * mu)

    n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))

    r0 = a * (1 - e * e) / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
    fact1 = n0 * math.tan(phi1) / r0

    _a1 = 500000 - easting
    dd0 = _a1 / (n0 * k0)
    fact2 = dd0 * dd0 / 2

    t0 = math.pow(math.tan(phi1), 2)
    Q0 = e1sq * math.pow(math.cos(phi1), 2)
    fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 - 9 * e1sq) * math.pow(dd0, 4) / 24

    fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 * e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720

    lof1 = _a1 / (n0 * k0)
    lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
    lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 * e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
    _a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
    _a3 = _a2 * 180 / math.pi

    latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi

    if not northernHemisphere:
        latitude = -latitude

    longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3

    return latitude, longitude


class ugrid:
    """
    A class for dealing with unstructured grid model output and converting to GNOME format
    Although I use variable names consistent with FVCOM, by passing in a var_map dict the 
    variable names can be customized for SELFE or ADCIRC
    
    Right now the attribute specifying whether the elements are orderd clockwise or counter
    clockwise needs to be manually added before writing to GNOME format (GNOME requres this, 
    but its not often specified in the model output)
        
    """
    
    def __init__(self,FileName=None):
            
        if FileName is not None:
            self.FileName = FileName
            if isinstance(FileName,list):
                self.Dataset = MFDataset(FileName)
            else:
                self.Dataset = Dataset(FileName)
            self.data = dict()
            self.atts = dict()
            
    def update(self,FileName):
        #point to a new nc file or url without reinitializing everything
        self.Dataset = Dataset(FileName)
            
    def get_dimensions(self,var_map):
        
        try:
            lat = self.Dataset.variables[var_map['latitude']]
            self.atts['lat'] = dict()
            for an_att in lat.ncattrs():
                self.atts['lat'][an_att] = getattr(lat,an_att)
            self.data['lat'] = lat[:]
            
            lon = self.Dataset.variables[var_map['longitude']]
            self.atts['lon'] = dict()
            for an_att in lon.ncattrs():
                self.atts['lon'][an_att] = getattr(lon,an_att)
            lon = lon[:]
            self.data['lon'] = (lon > 180).choose(lon,lon-360)
        except KeyError:
            print 'Lat/lon variables missing or named differently'
            pass
        
        try:
            time = self.Dataset.variables[var_map['time']]
            self.atts['time'] = dict()
            for an_att in time.ncattrs():
                self.atts['time'][an_att] = getattr(time,an_att)
            self.data['time'] = time[:]
            self.data['dtime'] = num2date(self.data['time'],time.units)
        except KeyError:
            pass
        
    def get_grid_topo(self,var_map):
        
        nv = self.Dataset.variables[var_map['nodes_surrounding_ele']]
        self.atts['nv'] = dict()
        for an_att in nv.ncattrs():
            self.atts['nv'][an_att] = getattr(nv,an_att)
        nv = nv[:]
#        if nv.min() == 0:
#            nv = nv + 1
        if nv.shape[0] > nv.shape[1]: #nv should be nv[num_faces,num_eles]
            self.data['nv'] = nv.transpose()
        else:
            self.data['nv'] = nv
            
        try: #should there be a check here for nbe starting at 0?
            nbe = self.Dataset.variables[var_map['eles_surrounding_ele']]
            self.atts['nbe'] = dict()
            for an_att in nbe.ncattrs():
                self.atts['nbe'][an_att] = getattr(nbe,an_att)               
            if nbe.shape[0] > nbe.shape[1]: #nbe should be nbe[num_faces,num_eles]
                self.data['nbe'] = nbe[:].transpose()
            else:
                self.data['nbe'] = nbe[:]
        except KeyError:
            print 'Building face-face connectivity'
            self.build_face_face_connectivity()
            
        try:
            #this returns ALL the edges (boundary and interior)
            edges = self.Dataset.variables[var_map['edge_node_connectivity']]
            self.atts['edges'] = dict()
            for an_att in edges.ncattrs():
                self.atts['edges'][an_att] = getattr(edges,an_att)  
            edges = edges[:]
#            if edges.min() == 0:
#                edges = edges + 1
            self.data['edges'] = edges
        except KeyError:
            print 'No edge information'
            #TODO: call code to determine edges here
            pass
        
    def get_data(self,var_map,tindex=None,nindex=None,zindex=0):
    
        ''' 
        var_map is a dict mapping model variable names to common names
        tindex can be used to subset in time --> tindex = [start,stop,step]
        nindex is for subsetting grid -- list of nodes/elements to get u/v on
        zindex is for z layer (FVCOM has surface at zindex = 0, SELFE zindex = -1)
        '''
    
        if tindex:
            self.data['time_ss'] = self.data['time'][tindex[0]:tindex[1]:tindex[2]]
        else:
            tindex = [0,len(self.data['time']),1]
                        
        u = self.Dataset.variables[var_map['u_velocity']]
        self.atts['u'] = dict()
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att)
        v = self.Dataset.variables[var_map['v_velocity']]
        self.atts['v'] = dict()
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att)
        
        if nindex is None:
            if len(u.shape)==3:
                self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],zindex,:]
                self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],zindex,:]    
            elif len(u.shape)==2:
                self.data['u'] = u[tindex[0]:tindex[1]:tindex[2],:]
                self.data['v'] = v[tindex[0]:tindex[1]:tindex[2],:]
            else:
                print "Error:velocity is not 2 or 3 dimensional"
                raise
        else: #Spatial subset -- under development but *mostly* working
            if u.shape[-1] == max(self.data['nbe'].shape):
                # velocities on elements
                id = np.where(np.diff(self.eles_in_ss) > 1)[0]
            else:
                # velocities on nodes
                id = np.where(np.diff(self.nodes_in_ss) > 1)[0]
            id2 = [-1]; id2.extend(id)

            #print 'Number of contiguous segments: ', len(id2)+1 
            firsttime = True
            for ii in range(len(id2)):
                #if not np.mod(ii,20):
                    #print ii,'of ', len(id2)
                sid = id2[ii]+1
                try:
                    fid = id2[ii+1]
                except IndexError:
                    fid = -1
                if  len(u.shape) == 3:
                    this_u = u[tindex[0]:tindex[1]:tindex[2],zindex,self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                    this_v = v[tindex[0]:tindex[1]:tindex[2],zindex,self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                elif len(u.shape)==2:
                    this_u = u[tindex[0]:tindex[1]:tindex[2],self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                    this_v = v[tindex[0]:tindex[1]:tindex[2],self.eles_in_ss[sid]-1:self.eles_in_ss[fid]]
                else:
                    print "Error:velocity is not 2 or 3 dimensional"
                    raise
                    
                if firsttime:
                    self.data['u'] = this_u.copy()
                    self.data['v'] = this_v.copy()
                    firsttime = False
                    ax = this_u.ndim-1
                else:  
                    self.data['u'] = np.concatenate((self.data['u'],this_u),axis=ax)
                    self.data['v'] = np.concatenate((self.data['v'],this_v),axis=ax)
        
        # sometimes these are numpy masked arrays -- GNOME can't deal
        if type(self.data['u']) is np.ma.core.MaskedArray:
            self.data['u'] = self.data['u'].data
        if type(self.data['v']) is np.ma.core.MaskedArray:
            self.data['v'] = self.data['v'].data
        
        #self.atts['v']['fill_value'] = self.atts['v']['missing_value']         
        #self.atts['u']['fill_value'] = self.atts['u']['missing_value']       
        
    def find_bndry_segs(self,subset=False):
    
        if subset:
            nv = self.data['nv_ss']
            nbe = self.data['nbe_ss']
        else:
            nv = self.data['nv']
            nbe = self.data['nbe']
            
        bnd = []
        for i, face in enumerate(nbe.transpose()):
            for j, neighbor in enumerate(face):
                if neighbor == 0:
                    if j == 0:
                        bound = [nv[1,i], nv[2,i]]
                    elif j == 1:
                        bound = [nv[2,i], nv[0,i]]
                    elif j == 2:
                        bound = [nv[0,i], nv[1,i]]
#                    if j == num_vertices-1:
#                        bound = [nv[-1,i], nv[0,i]]
#                    else:
#                        bound = [nv[j,i], nv[j+1,i]]
                    bnd.append(bound)
        
        return bnd

    def order_boundary(self,b,seg_types):
    
        obnd = dict()
        otype = dict()
        bnd_number = 0
        obnd[bnd_number] = [b.pop(0),]
        otype[bnd_number] = [seg_types.pop(0),]
        while len(b)>0:
            idx = [i for i,edge in enumerate(b) if edge[0]==obnd[bnd_number][-1][-1]]
            if len(idx) == 1:
                obnd[bnd_number].append(b.pop(idx[0]))
                otype[bnd_number].append(seg_types.pop(idx[0]))
            else:
                bnd_number = bnd_number + 1
                obnd[bnd_number] = [b.pop(0),]
                otype[bnd_number] = [seg_types.pop(0),]
                
        #format for GNOME ([node1,node2,bnd_num,bnd_type] - bnd_type=1 for open, 2 for closed)
        boundary = []
        for i, a_bnd in obnd.iteritems():
            for j, seg in enumerate(a_bnd):
                #TODO -- need to make separate method for adding 1 to nv,boundary...
                boundary.append([seg[0], seg[1], i, otype[i][j]])
            
        self.data['bnd'] = boundary
        self.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'}  
    
    def write_bndry_file(self,bnd_file):
        pass
    
    def read_bndry_file(self,bnd_file): 
 
        bnd = []
        f = open(bnd_file,'r')
        for line in f:
            vals = [int(val) for val in line.split()]
            bnd.append(vals)
        
        self.data['bnd'] = np.array(bnd)
        self.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'}    
        
    def write_unstruc_grid(self,ofn):
        
        """
        
        Write GNOME compatible netCDF file (netCDF3) from unstructured (triangular) grid data
        
        """  
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','Triangular')
        
        # test u/v dimensions
        if self.data['u'].shape != self.data['v'].shape:
            print 'u/v dimensions differ'
            raise
        
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError:
            if self.data['u'].shape[0] != len(self.data['time']):
                print 'Dimensions of u/v do not match time variable'
                raise
                
        lon_key = 'lon'; lat_key = 'lat'
        nv_key = 'nv'; nbe_key = 'nbe'
        # determine if its a subset of the grid
        try:
            if self.data['u'].shape[-1] == len(self.data['lon_ss']) or \
                self.data['u'].shape[-1] == self.data['nbe_ss'].shape[-1]:
                lon_key = 'lon_ss'; lat_key = 'lat_ss'
                nv_key = 'nv_ss'; nbe_key = 'nbe_ss'
        except KeyError:
            if self.data['u'].shape[-1] != len(self.data['lon']) and \
                self.data['u'].shape[-1] != self.data['nbe'].shape[-1]:
                print 'Dimensions of u/v do not match grid variables'
                raise 
                
        # add Dimensions
        nc.createDimension('time',None)
        nc.createDimension('node',len(self.data[lon_key]))
        nc.createDimension('nele',np.shape(self.data[nbe_key])[1])
        nc.createDimension('nbnd',len(self.data['bnd']))
        nc.createDimension('nbi',4)
        nc.createDimension('three',3)
        #nc.createDimension('sigma',1) #coming soon?
        
        try:
            ufill = self.atts['u']['_FillValue']
            vfill = self.atts['v']['_FillValue']
        except KeyError:
            try:
                ufill = self.atts['u']['missing_value']
                vfill = self.atts['v']['missing_value']
            except KeyError:
                ufill = 999.
                vfill = 999.
        
        # create variables
        nc_time = nc.createVariable('time','f4',('time',))
        nc_lon = nc.createVariable('lon','f4',('node'))
        nc_lat = nc.createVariable('lat','f4',('node'))
        nc_nbe = nc.createVariable('nbe','int32',('three','nele'))
        nc_nv = nc.createVariable('nv','int32',('three','nele'))
        nc_bnd = nc.createVariable('bnd','int32',('nbnd','nbi'))
        
        if self.data['u'].shape[-1] == len(self.data[lon_key]): #velocities on nodes
            nc_u = nc.createVariable('u','f4',('time','node'),fill_value=ufill)
            nc_v = nc.createVariable('v','f4',('time','node'),fill_value=vfill)
        else: #velocities on elements
            nc_u = nc.createVariable('u','f4',('time','nele'),fill_value=ufill)
            nc_v = nc.createVariable('v','f4',('time','nele'),fill_value=vfill)
        
        #adjust time if necessary
        ref_time = self.atts['time']['units'].split(' since ')[1]
        ref_year = int(ref_time[0:4])
        if ref_year < 1970:
            print 'Adjusting reference time'
            self.data[t_key],self.atts['time']['units'] = \
                nctools.adjust_time(self.data[t_key],self.atts['time']['units'])
        
        #add data to netcdf file
        nc_time[:] = self.data[t_key]
        nc_lon[:] = self.data[lon_key]
        nc_lat[:] = self.data[lat_key]
        nc_u[:] = self.data['u']
        nc_v[:] = self.data['v']
        nc_bnd[:] = self.data['bnd']
        nc_nbe[:] = self.data[nbe_key]
        nc_nv[:] = self.data[nv_key]
        
        #add variable attributes to netcdf file
        for an_att in self.atts['time'].iteritems():
           setattr(nc_time,an_att[0],an_att[1])
        
        for an_att in self.atts['lon'].iteritems():
            setattr(nc_lon,an_att[0],an_att[1])
        
        for an_att in self.atts['lat'].iteritems():
            setattr(nc_lat,an_att[0],an_att[1])
        
        for an_att in self.atts['bnd'].iteritems():
            setattr(nc_bnd,an_att[0],an_att[1])

        for an_att in self.atts['nbe'].iteritems():
            setattr(nc_nbe,an_att[0],an_att[1])

        for an_att in self.atts['nv'].iteritems():
            setattr(nc_nv,an_att[0],an_att[1])
        
        for an_att in self.atts['u'].iteritems():
            if an_att[0] != '_FillValue':
                setattr(nc_u,an_att[0],an_att[1])
        
        for an_att in self.atts['v'].iteritems():
            if an_att[0] != '_FillValue':
                setattr(nc_v,an_att[0],an_att[1])
        
        nc.close()
    
    def write_unstruc_grid_only(self,ofn):
        
        """
        
        Write netCDF file (netCDF3) of grid variables for unstructured (triangular) grid
        
        """  
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','Triangular')

        # add Dimensions
        nc.createDimension('node',len(self.data['lon']))
        nc.createDimension('nele',np.shape(self.data['nbe'])[1])
        nc.createDimension('nbnd',len(self.data['bnd']))
        nc.createDimension('nbi',4)
        nc.createDimension('three',3)
        #nc.createDimension('sigma',1) #coming soon?
        
        # create variables
        nc_lon = nc.createVariable('lon','f4',('node'))
        nc_lat = nc.createVariable('lat','f4',('node'))
        nc_nbe = nc.createVariable('nbe','int32',('three','nele'))
        nc_nv = nc.createVariable('nv','int32',('three','nele'))
        nc_bnd = nc.createVariable('bnd','int32',('nbnd','nbi'))
        
        #add data to netcdf file
        nc_lon[:] = self.data['lon']
        nc_lat[:] = self.data['lat']
        nc_bnd[:] = self.data['bnd']
        nc_nbe[:] = self.data['nbe']
        nc_nv[:] = self.data['nv']
        
        for an_att in self.atts['lon'].iteritems():
            setattr(nc_lon,an_att[0],an_att[1])   
        for an_att in self.atts['lat'].iteritems():
            setattr(nc_lat,an_att[0],an_att[1])       
        for an_att in self.atts['bnd'].iteritems():
            setattr(nc_bnd,an_att[0],an_att[1])
        for an_att in self.atts['nbe'].iteritems():
            setattr(nc_nbe,an_att[0],an_att[1])
        for an_att in self.atts['nv'].iteritems():
            setattr(nc_nv,an_att[0],an_att[1])
        
        nc.close()
        
    def find_nodes_eles_in_ss(self,nl,sl,wl,el):
       
        print 'Total number of eles: ', self.data['nbe'].shape[1]
        print 'Total number of nodes: ', self.data['lon'].shape[0]
        
        #returns lists of eles and nodes, plus truncated and edited topology arrays (nbe, nv)
        subset_lat = np.nonzero(np.logical_and(self.data['lat']>=sl,self.data['lat']<=nl))[0]
        subset_lon = np.nonzero(np.logical_and(self.data['lon']>=wl,self.data['lon']<=el))[0]
        self.nodes_in_ss = np.intersect1d(subset_lat,subset_lon) + 1 #node numbering starts at 1 (subtract one for indexing lon/lat)
        
        self.eles_in_ss = []
        nv_ss = []; nbe_ss = []; 
        
        #determine which nodes are in subset boundary and elements with all nodes in ss
        print 'Finding nodes and entire elements in ss'
        for ii, ele in enumerate(self.data['nv'].transpose()):
            #if all of the nodes are in subset domain keep -- otherwise get rid of it
            if (self.nodes_in_ss == ele[0]).any() and (self.nodes_in_ss == ele[1]).any() \
             and (self.nodes_in_ss == ele[2]).any():
                nv_ss.append(ele)
                nbe_ss.append(self.data['nbe'][:,ii])
                self.eles_in_ss.append(ii+1) #ele numbering starts at 1
            else:
                pass
            
        print 'Number of eles in ss: ', len(self.eles_in_ss)
        print 'Number of nodes in ss: ', len(self.nodes_in_ss)
              
        nbe_ss = np.array(nbe_ss).transpose()
        nv_ss = np.array(nv_ss).transpose()
        self.eles_in_ss = np.array(self.eles_in_ss)

        print 'Remapping nodes and elements'
        #now remap nbe_ss, nv_ss to number of remaining nodes, and elements
        nv_ssr = nv_ss.copy()
        nbe_ssr = nbe_ss.copy()
        for ii in range(len(self.eles_in_ss)):
            for jj in range(3):
                nid = np.searchsorted(self.nodes_in_ss, nv_ss[jj,ii], side='left')
                nv_ssr[jj,ii] = nid+1
                if nbe_ss[jj,ii] != 0:
                    eid = np.searchsorted(self.eles_in_ss, nbe_ss[jj,ii], side='left')
                    if eid >= len(self.eles_in_ss) or self.eles_in_ss[eid] != nbe_ss[jj,ii]:
                        nbe_ssr[jj,ii] = 0
                    else:
                        nbe_ssr[jj,ii] = eid+1
                                   
        self.data['nbe_ss'] = nbe_ssr
        self.data['nv_ss'] = nv_ssr
        self.data['lon_ss'] = self.data['lon'][self.nodes_in_ss-1]
        self.data['lat_ss'] = self.data['lat'][self.nodes_in_ss-1]
        
    def remap_bry_nodes(self,bndry_file):
      
        #find all the land segments and re-number to match new subset node numbers
        #The outer boundary will be determined in write_bndry file then this segment
        #info will be used to mark land segments
        
        print 'Remapping boundary segs to new subset node numbers'
        f = open(bndry_file)
        self.ss_land_bry_segs = []
        for line in f:
            node1,node2,bnumber,flag = map(int,line.split())
            node1_id = np.where(self.nodes_in_ss == node1)[0]
            node2_id = np.where(self.nodes_in_ss == node2)[0]
            if len(node1_id) > 0 and len(node2_id) > 0 and flag == 0:
                self.ss_land_bry_segs.append([node1_id[0]+1,node2_id[0]+1,flag])
        f.close()
    
    
    def build_face_face_connectivity(self):
        """
        builds the triangular connectivity array (nbe)
        essentially giving the neighbors of each triangle (face)
        """       
                
        num_vertices = 3
        num_faces = max(self.data['nv'].shape)
        face_face = np.zeros( (num_faces, num_vertices), dtype=np.int32  )
        face_face += -1 # fill with -1

        # loop through all the triangles to find the matching edges:
        edges = {} # dict to store the edges in 
        for i, face in enumerate(self.data['nv'].transpose()):
            # loop through edges:
            for j in range(num_vertices):
                edge = (face[j-1], face[j])
                if edge[0] > edge[1]: # sort the node numbers
                    edge = (edge[1], edge[0]) 
                # see if it is already in there
                prev_edge = edges.pop(edge, None)
                if prev_edge is not None:
                    face_num, edge_num = prev_edge
                    face_face[i,j] = face_num
                    face_face[face_num, edge_num] = i
                else:
                    edges[edge] = (i, j)
        nbe = (face_face + 1)
        nbe = nbe[:,[2,0,1]].transpose()
        self.data['nbe'] = nbe
        self.atts['nbe'] = {'long_name': 'elements surrounding element'}

class cgrid():
    
    def __init__(self,FileName=None):
        
        if FileName is not None:
            self.FileName = FileName
            if isinstance(FileName,list):
                self.Dataset = MFDataset(FileName)
            else:
                self.Dataset = Dataset(FileName)
            self.data = dict()
            self.atts = dict()
            self.grid = dict()
            
    def get_dimensions(self,var_map):
        
        self.time = self.Dataset.variables[var_map['time']]  
        self.atts['time'] = {}
        for an_att in self.time.ncattrs():
            self.atts['time'][an_att] = getattr(self.time,an_att) 
        self.data['time'] = self.time[:]
    
        lon = self.Dataset.variables[var_map['lon']]
        self.atts['lon'] = {}
        for an_att in lon.ncattrs():
            self.atts['lon'][an_att] = getattr(lon,an_att)
        self.data['lon'] = lon[:]
        
        lat = self.Dataset.variables[var_map['lat']]
        self.atts['lat'] = {}
        for an_att in lat.ncattrs():
            self.atts['lat'][an_att] = getattr(lat,an_att)
        self.data['lat'] = lat[:]        
    
    def subset_pt_in_poly(self,bbox,stride=1,lat='lat',lon='lon'):
        '''
        bbox = [slat,wlon,nlat,elon]
        can pass in lat/lon names to specify which grid the subset is done on (for c-grids)
        '''
        glat = self.data[lat]
        glon = self.data[lon]
        
        bbox=np.array(bbox)
        mypath=np.array([bbox[[1,3,3,1]],bbox[[0,0,2,2]]]).T
        p = path.Path(mypath)
        points = np.vstack((glon.flatten(),glat.flatten())).T   
        n,m = np.shape(glon)
        inside = p.contains_points(points).reshape((n,m))
        ii,jj = np.meshgrid(xrange(m),xrange(n))
        self.x = [min(ii[inside]),max(ii[inside])+1,stride]
        self.y = [min(jj[inside]),max(jj[inside])+1,stride]
        
    def subset(self,bbox,stride=1,dl=0,lat='lat',lon='lon'):
        '''
        bbox = [slat,wlon,nlat,elon]
        can pass in lat/lon names to specify which grid the subset is done on (for c-grids)
        '''
        glat = self.data[lat]
        glon = self.data[lon]
        
        sl = bbox[0]
        nl = bbox[2]
        wl = bbox[1]
        el = bbox[3]
        
        if (abs(np.nanmax(glat)-nl) < 1e-3) and (abs(np.nanmin(glon)-wl) < 1e-3): #original values
            self.y = [0,np.size(glat,0),1]
            self.x = [[0,np.size(glat,1),1]]
        else: #do subset
            
            if dl == 0:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_and(glon>=wl,glon<=el)))
            else:
                [yvec,xvec] = np.where(np.logical_and(np.logical_and(glat>=sl,glat<=nl),np.logical_or(glon>=wl,glon<=el)))
            
            if len(yvec) > 2 and len(xvec) > 2:
                y1 = min(yvec)
                y2 = max(yvec)+1
                x1 = min(xvec)
                x2 = max(xvec)+1
                self.y = [y1,y2,stride]
                self.x = [x1,x2,stride]           
            else:
                self.y = [0,np.size(glat,0),1]
                self.x = [[0,np.size(glat,1),1]]
               
    def get_grid_info(self,grid_vars=['mask'],yindex=None,xindex=None):
    
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
        
        for var in grid_vars: 
        #these are 2D spatial vars (need to add vertical info) for 3d
            try:
                self.grid[var] = self.Dataset.variables[var][y1:y2:step,x1:x2:step]
            except KeyError:
                pass

    
    def get_data(self,var_map,tindex=None,yindex=None,xindex=None,zindex=0,is3d=False):
    
        ''' 
        var_map is a dict mapping model variable names to common names
        tindex can be used to subset in time --> tindex = [start,stop,step]
        xindex and yindex are for subsetting grid 
        zindex is for z layer (surface at zindex = 0 or -1)
        
        '''
        if tindex is None:
            self.data['time_ss'] = self.data['time']
            t1 = 0; t2 = len(self.data['time']); ts = 1
        else:
            t1 = tindex[0]; t2 = tindex[1]; ts = tindex[2]
            self.data['time_ss'] = self.data['time'][t1:t2:ts]
            
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
            self.data['lon_ss'] = self.data['lon'][y1:y2:step,x1:x2:step]
            self.data['lat_ss'] = self.data['lat'][y1:y2:step,x1:x2:step]
        
        u = self.Dataset.variables[var_map['u']]
        self.atts['u'] = {}
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att) 

        v = self.Dataset.variables[var_map['v']]
        self.atts['v'] = {}
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att) 
        
        self.data['u'] = u[t1:t2:ts,zindex,y1:y2:step,x1:x2:step]
        self.data['v'] = v[t1:t2:ts,zindex,y1:y2:step,x1:x2:step]

        
    def write_nc(self,ofn,is3d=False):
        """
        Write GNOME compatible netCDF file (netCDF3)
        * velocities are on center points and rotated to north/east
        * lat/lon can EITHER have the same dimensions as u/v or be one larger
          in both x/y dimensions (i.e. it defines the stencil) - GNOME will
          treat these two cases differently
        
        """
        
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','curvilinear')
    
        # test u/v dimensions
        if self.data['u'].shape != self.data['v'].shape:
            raise Exception('u/v dimensions differ')
        
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError: #TODO -- if it has key time_ss but it doesn't match that case is not caught
            if self.data['u'].shape[0] != len(self.data['time']):
                raise Exception('Dimensions of u/v do not match time variable')
                
        lon_key = 'lon'; lat_key = 'lat'
        # determine if its a subset of the grid
        try:
            lon_ss_shape = self.data['lon_ss'].shape
            lon_ss_shape_red = (lon_ss_shape[0]-1,lon_ss_shape[1]-1)
            if self.data['u'].shape[-2:] == lon_ss_shape or \
                self.data['u'].shape[-2:] == lon_ss_shape_red:
                lon_key = 'lon_ss'; lat_key = 'lat_ss'
        except KeyError:
            lon_shape = self.data['lon'].shape
            lon_shape_red = (lon_shape[0]-1,lon_shape[1]-1)
            if self.data['u'].shape[-2:] != lon_shape and \
                self.data['u'].shape[-2:] != lon_shape_red:
                raise Exception('Dimensions of u/v do not match grid variables')
        
	                
        x = self.data[lon_key].shape[1]
        y = self.data[lat_key].shape[0]
        
        if self.data[lon_key].shape == self.data['u'].shape[-2:]:
            xc = x
            yc = y
        else:  
            xc = x-1
            yc = y-1

        # add Dimensions
        nc.createDimension('x',x)
        nc.createDimension('y',y)
        nc.createDimension('xc',xc)
        nc.createDimension('yc',yc)
        nc.createDimension('time',None)
    
        try:
            ufill = self.atts['u']['_FillValue']
            vfill = self.atts['v']['_FillValue']
        except KeyError:
            try:
                ufill = self.atts['u']['missing_value']
                vfill = self.atts['v']['missing_value']
            except KeyError:
                ufill = 999.
                vfill = 999.
    
        # create variables
        nc_time = nc.createVariable('time','f4',('time',))
        nc_lon = nc.createVariable('lon','f4',('y','x'))
        nc_lat = nc.createVariable('lat','f4',('y','x'))
        
        if self.atts.has_key('wind'):
            nc_u = nc.createVariable('air_u','f4',('time','yc','xc'), \
                fill_value=ufill)
            nc_v = nc.createVariable('air_v','f4',('time','yc','xc'), \
                fill_value=vfill)
        elif is3d:
            pass
        else:
            nc_u = nc.createVariable('water_u','f4',('time','yc','xc'), \
                fill_value=ufill)
            nc_v = nc.createVariable('water_v','f4',('time','yc','xc'), \
                fill_value=vfill)
        
         # add data
        nc_lon[:] = self.data[lon_key]
        nc_lat[:] = self.data[lat_key]
        
        #t_key='time'
        #!!!!!!!!!!!!Add 3d
        if len(self.data['u'].shape) == 2:
            nc_time[0] = self.data[t_key]
            nc_u[0,:] = self.data['u']
            nc_v[0,:] = self.data['v']
        else:
            nc_time[:] = self.data[t_key]
            nc_u[:] = self.data['u']
            nc_v[:] = self.data['v']
            
        if self.grid.has_key('mask'):
            nc_mask = nc.createVariable('mask','f4',('yc','xc'))
            nc_mask[:] = self.grid['mask']
        
        #pdb.set_trace()
        # add variable attributes from 'atts' (nested dict object)
        for key,val in self.atts['time'].iteritems():
            if not key.startswith('_'):
                setattr(nc_time,key,val)
            
        for key,val in self.atts['u'].iteritems():
            if not key.startswith('_'):
                setattr(nc_u,key,val)
    
        for an_att in self.atts['v'].iteritems():
            if not key.startswith('_'):
                setattr(nc_v,key,val)
    
        nc.close()
    
class roms(cgrid):
    """
    A class for dealing with curvilinear grid model output and converting to GNOME format
    Requires passing in a var_map dict so the 
    variable names can be customized for different models or datasets
            
   """
         
    def get_dimensions(self,var_map):
        
        #import pdb
        #nc=Dataset('http://barataria.tamu.edu:8080/thredds/dodsC/txla/roms/2014/ocean_his_03.nc','r')
        #nc=Dataset('txla_subset_HIS_Mar2014.nc','r')
        #print nc.variables
        #print nc.variables['lat_rho']
        #pdb.set_trace()
        self.time = self.Dataset.variables[var_map['time']]
        self.atts['time'] = {}
        for an_att in self.time.ncattrs():
            self.atts['time'][an_att] = getattr(self.time,an_att) 
        self.data['time'] = self.time[:]
        
        #load lat/lon for rho, u, and v grids
        for var in ['lat_rho','lat_u','lat_v','lon_rho','lon_u','lon_v']:
            ds_var = self.Dataset.variables[var]
            self.atts[var] = {}
            for an_att in ds_var.ncattrs():
                self.atts[var][an_att] = getattr(ds_var,an_att)
            self.data[var] = ds_var[:]
        
        #Now load or create P grid lat/lon (sometimes not included in ROMS output)
        try:
            lon_psi = self.Dataset.variables['lon_psi']
            self.atts['lon_psi'] = {}            
            for an_att in lon_psi.ncattrs():
                self.atts['lon_psi'][an_att] = getattr(lon_psi,an_att) 
            self.data['lon_psi'] = lon_psi[:]
            lat_psi = self.Dataset.variables['lat_psi']
            self.atts['lat_psi'] = {}
            for an_att in lat_psi.ncattrs():
                self.atts['lat_psi'][an_att] = getattr(lat_psi,an_att) 
            self.data['lat_psi'] = lat_psi[:]
        except KeyError:
            self.data['lon_psi'] = (self.data['lon_rho'][0:-1,0:-1]+self.data['lon_rho'][1:,1:])*0.5
            self.atts['lon_psi'] = self.atts['lon_rho']
            self.atts['lon_psi']['long_name'] = 'longitude of PSI-points'
            self.data['lat_psi'] = (self.data['lat_rho'][0:-1,0:-1]+self.data['lat_rho'][1:,1:])*0.5
            self.atts['lat_psi'] = self.atts['lat_rho']
            self.atts['lat_psi']['long_name'] = 'latitude of PSI-points'
                

    def get_grid_info(self,yindex=None,xindex=None,is3d=False):
        
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon_psi'].shape[1]
            y1 = 0; y2 = self.data['lon_psi'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]
            x1 = xindex[0]; x2 = xindex[1]
        
	
        self.grid['mask'] = self.Dataset.variables['mask_rho'][y1:y2+1,x1:x2+1] 
        self.grid['angle'] = self.Dataset.variables['angle'][y1:y2+1,x1:x2+1] 
        
        if is3d:
            self.grid['h'] = self.Dataset.variables['h'][y1:y2+1,x1:x2+1] 
            self.grid['hc'] = self.Dataset.variables['hc'][:]
            self.grid['Cs_r'] = self.Dataset.variables['Cs_r'][:]
            self.grid['sc_r'] = self.Dataset.variables['s_rho'][:]
    
    def get_data(self,var_map,tindex=None,yindex=None,xindex=None,is3d=False,interp=True):
        
        '''
        In this case, lon/lat on psi (P) grid, u on u-grid, v on v-grid
        
        '''
        if tindex is None:
            self.data['time_ss'] = self.data['time']
            t1 = 0; t2 = len(self.data['time']); ts = 1
        else:
            t1 = tindex[0]; t2 = tindex[1]; ts = tindex[2]
            self.data['time_ss'] = self.data['time'][t1:t2:ts]
            
        if xindex is None and yindex is None:
            x1 = 0; x2 = self.data['lon_psi'].shape[1]; step = 1
            y1 = 0; y2 = self.data['lon_psi'].shape[0]
        else:
            y1 = yindex[0]; y2 = yindex[1]; step = yindex[2]
            x1 = xindex[0]; x2 = xindex[1]
	    if y2==self.data['lon_psi'].shape[0] or x2==self.data['lon_psi'].shape[1]:
            	## The larger index reach the maximum
        	self.data['lon_psi_ss'] = self.data['lon_psi'][y1:y2+1:step,x1:x2:step]
                self.data['lat_psi_ss'] = self.data['lat_psi'][y1:y2+1:step,x1:x2:step]
	    else:
            	self.data['lon_psi_ss'] = self.data['lon_psi'][y1:y2+1:step,x1:x2+1:step]
            	self.data['lat_psi_ss'] = self.data['lat_psi'][y1:y2+1:step,x1:x2+1:step]
        
        u = self.Dataset.variables['u']
        self.atts['u'] = {}
        for an_att in u.ncattrs():
            self.atts['u'][an_att] = getattr(u,an_att) 
        v = self.Dataset.variables['v']
        self.atts['v'] = {}
        for an_att in v.ncattrs():
            self.atts['v'][an_att] = getattr(v,an_att) 
	
        if is3d: 
            u_on_upts = u[t1:t2+1:ts,:,y1:y2+1,x1:x2]
            v_on_vpts = v[t1:t2+1:ts,:,y1:y2,x1:x2+1]
        else:
	    if u.ndim == 4:
                u_on_upts = u[t1:t2+1:ts,-1,y1:y2+1,x1:x2]
                v_on_vpts = v[t1:t2+1:ts,-1,y1:y2,x1:x2+1]
	    else:
		u_on_upts = u[t1:t2+1:ts,y1:y2+1,x1:x2]
                v_on_vpts = v[t1:t2+1:ts,y1:y2,x1:x2+1]

        #replace nans or fill values with 0 for interpolating to rho grid
        u_on_upts = (np.isnan(u_on_upts)).choose(u_on_upts,0)
        v_on_vpts = (np.isnan(v_on_vpts)).choose(v_on_vpts,0)
        u_on_upts = (u_on_upts > 1e10).choose(u_on_upts,0)
        v_on_vpts = (v_on_vpts > 1e10).choose(v_on_vpts,0)
               
        if interp:
            self.data['u'],self.data['v'] = self.interp_and_rotate(u_on_upts,v_on_vpts)
        else:
            self.data['u'] = u_on_upts
            self.data['v'] = v_on_vpts

            
    
    def interp_and_rotate(self,u,v,is3d=False):
        '''Calculate u/v on rho points -- we lose exterior most u/v values
        Then rotate to north/east
        '''
        self.grid['angle'] = self.grid['angle'][1:-1,1:-1]
        self.grid['mask'] = self.grid['mask'][1:-1,1:-1]
        cosa = (np.cos(self.grid['angle']) * self.grid['mask'])
        sina = (np.sin(self.grid['angle']) * self.grid['mask'])
        if is3d:
            u_rot = np.zeros([u.shape[0],u.shape[1],v.shape[2]-1,v.shape[2]-1])
            v_rot = np.zeros_like(u_rot)
            for z in u.shape[1]:
                u_rho = (u[:,z,1:-1,:-1] +u[:,z,:-1,1:])/2. 
                v_rho = (v[:,z,:-1,1:-1] + v[:,z,1:,1:-1])/2.
                u_rot[:,z,:,:] = u_rho * cosa - v_rho * sina
                v_rot[:,z,:,:] = u_rho * sina + v_rho * cosa
        else:
            u_rho = (u[...,1:-1,:-1] + u[...,1:-1,1:])/2. 
            #u_rho = (u[...,1:,:-1] + u[...,1:,1:])/2. 
            v_rho = (v[...,:-1,1:-1] + v[...,1:,1:-1])/2.
            #v_rho = (v[...,:-1,1:] + v[...,1:,1:])/2.
            #pdb.set_trace()
            #u_rho=u_rho[0,:,:]
            #v_rho=v_rho[0,:,:]
            #cosa=cosa[0:-1,:]
            #sina=sina[0:-1,:]
            u_rot = u_rho * cosa - v_rho * sina
            v_rot = u_rho * sina + v_rho * cosa
            
        return u_rot, v_rot
        
    def reduce_latlon_mesh_for_GNOME(self):
              
        if self.data.has_key('lon_psi_ss'): #subset
            self.data['lon_psi_ss'] = self.data['lon_psi_ss'][:-1,:-1]
            self.data['lat_psi_ss'] = self.data['lat_psi_ss'][:-1,:-1]   
        else:
            self.data['lon_psi'] = self.data['lon_psi'][:-1,:-1]
            self.data['lat_psi'] = self.data['lat_psi'][:-1,:-1]

    def write_nc_native(self,ofn,is3d=False):
        """
      
        Write GNOME compatible netCDF file (netCDF3)
        Maintain u and v on sepearate grids
        
        %TODO: once this is implemented in GNOME, rename to "write_nc"
        to replace method in base class
        
        
        """
        nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
        
        # Global Attributes
        setattr(nc,'grid_type','curvilinear')
        setattr(nc,'grid_case','arakawa-c')
    
        x = self.data['lon_psi'].shape[1]
        y = self.data['lat_psi'].shape[0]

        xc = x+1 #rho on center points plus outside boundary
        yc = y+1
        
        # determine if its a subset in time
        t_key = 'time'
        try:
            if self.data['u'].shape[0] == len(self.data['time_ss']):
                t_key = 'time_ss'
        except KeyError:
            if self.data['u'].shape[0] != len(self.data['time']):
                raise Exception('Dimensions of u/v do not match time variable')
        
        
        
        
        
        # add Dimensions
        nc.createDimension('x',x)
        nc.createDimension('y',y)
        nc.createDimension('xc',xc)
        nc.createDimension('yc',yc)
        nc.createDimension('time',None)
    
        try:
            ufill = self.atts['u']['_FillValue']
            vfill = self.atts['v']['_FillValue']
        except KeyError:
            try:
                ufill = self.atts['u']['missing_value']
                vfill = self.atts['v']['missing_value']
            except KeyError:
                ufill = 999.
                vfill = 999.
    
        # create variables
        nc_time = nc.createVariable('time','f4',('time',))
        nc_lonp = nc.createVariable('lon_psi','f4',('y','x'))
        nc_latp = nc.createVariable('lat_psi','f4',('y','x'))
        nc_lonr = nc.createVariable('lon_rho','f4',('yc','xc'))
        nc_latr = nc.createVariable('lat_rho','f4',('yc','xc'))
        nc_lonu = nc.createVariable('lon_u','f4',('yc','x'))
        nc_latu = nc.createVariable('lat_u','f4',('yc','x'))
        nc_lonv = nc.createVariable('lon_v','f4',('y','xc'))
        nc_latv = nc.createVariable('lat_v','f4',('y','xc'))
        nc_angle = nc.createVariable('angle','f4',('yc','xc'))

        if is3d:
            pass
        else:
            nc_u = nc.createVariable('water_u','f4',('time','yc','x'), \
                fill_value=ufill)
            nc_v = nc.createVariable('water_v','f4',('time','y','xc'), \
                fill_value=vfill)
        
         # add data
        nc_lonp[:] = self.data['lon_psi']
        nc_latp[:] = self.data['lat_psi']
        nc_lonr[:] = self.data['lon_rho']
        nc_latr[:] = self.data['lat_rho']
        nc_lonu[:] = self.data['lon_u']
        nc_latu[:] = self.data['lat_u']
        nc_lonv[:] = self.data['lon_v']
        nc_latv[:] = self.data['lat_v']
        nc_angle[:] = self.grid['angle']
        
        #!!!!!!!!!!!!Add 3d
        if len(self.data['u'].shape) == 2:
            nc_time[0] = self.data[t_key]
            nc_u[0,:] = self.data['u']
            nc_v[0,:] = self.data['v']
        else:
            nc_time[:] = self.data[t_key]
            nc_u[:] = self.data['u']
            nc_v[:] = self.data['v']
            
        if self.grid.has_key('mask'):
            nc_mask = nc.createVariable('mask_rho','f4',('yc','xc'))
            nc_mask[:] = self.grid['mask']

        pdb.set_trace()
        # add variable attributes from 'atts' (nested dict object)
        for key,val in self.atts['time'].iteritems():
            if not key.startswith('_'):
                setattr(nc_time,key,val)
            
        for key,val in self.atts['u'].iteritems():
            if not key.startswith('_'):
                setattr(nc_u,key,val)
    
        for an_att in self.atts['v'].iteritems():
            if not key.startswith('_'):
                setattr(nc_v,key,val)
    
    
        nc.close()

