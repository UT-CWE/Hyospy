# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:55:05 2017

@author: dongyu
"""

from upper_wrapper import upper_wrapper

import pdb

starttime='2017-05-28-00'
endtime='2017-05-29-00'

UW = upper_wrapper()

## Hydrodynamic model
UW.hydro_model = 'ROMS'  # hydrodynamic model options: 'SUNTANS', 'ROMS', 'BLENDED'

## ROMS datasource
UW.ROMS_datasource = 'online' # ROMS data source options: 'online', 'offline'

## Oil spill model
UW.runGNOME = True
UW.runTracPy = False

## Probability Map
UW.probability_map = True
UW.google_earth = False
UW.mpl = 8 # probability calculation time interval

## ensemble parameters
UW.number = 1
UW.interval = 3600 # units: seconds

## SUNTANS IC and BC
UW.OBC_opt = 'file'  # Type 3 boundary condition option: 'constant',
#'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'
UW.IC_opt = 'constant'

#UW(starttime, endtime, 20, init_latlon=[29.463089, -94.843460]) #SUNTANS domain
UW(starttime, endtime, 20, init_latlon=[28.353786, -95.315109])  #ROMS domain



