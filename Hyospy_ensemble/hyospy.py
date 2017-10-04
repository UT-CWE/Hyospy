# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:55:05 2017

@author: dongyu
"""

from upper_wrapper import upper_wrapper

import pdb

## Hydrodynamic model time info
starttime='2017-05-29-00'
endtime='2017-05-31-00'

## Oil spill model time info
starttime2='2017-05-29-20'  # starttime for oil spill model
period=20

UW = upper_wrapper()

## Hydrodynamic model
UW.hydro_model = 'ROMS'  # hydrodynamic model options: 'SUNTANS', 'ROMS', 'BLENDED'
UW.hydro_run = False  # choose whether or not to run a hydrodynamic model first

## ROMS datasource
UW.ROMS_datasource = 'offline' # ROMS data source options: 'online', 'offline'

## Oil spill model
UW.runGNOME = True
UW.runTracPy = False

## GNOME settings
UW.gnome_subset = True       # For blended current product, "subset=True" makes GNOME run faster
UW.gnome_bbox = [28.17,-95.53,30.0,-93.9]

## Probability Map
UW.probability_map = False
UW.google_earth = False
UW.mpl = 8 # probability calculation time interval

## ensemble parameters
UW.number = 1
UW.interval = 500 # units: seconds

## SUNTANS IC and BC
UW.OBC_opt = 'ROMSFILE'  # Type 3 boundary condition option: 'constant',
#'file','OTIS', 'ROMS', 'ROMSOTIS','ROMSFILE', 'ROMSOTISFILE'
UW.IC_opt = 'ROMS'

#UW(starttime, endtime, starttime2, period, init_latlon=[29.463089, -94.843460]) #SUNTANS domain
UW(starttime, endtime, starttime2, period, init_latlon=[28.353786, -95.315109])  #ROMS domain



