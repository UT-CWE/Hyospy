# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:09:41 2017

@author: dongyu
"""

import lib
import lib.SUNTANS
from DownloadTool import downloadROMS

import pdb

#starttime = '2017-05-01-11'
#endtime = '2017-05-04-00'

starttime = '2017-05-03-00'
endtime = '2017-05-05-00'

downloadROMS(starttime, endtime, ROMSsrc='forecast')
