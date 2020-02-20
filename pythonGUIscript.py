# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 10:21:54 2020

@author: mitch
"""

import numpy as np
import astropy as ast
import eleanor
from astropy import units as u
from astropy.coordinates import SkyCoord
from tkinter import *
import gzip
from astropy.io import fits
from astropy.table import table
from astroquery.mast import Tesscut, Observations

#Code cannot be used, the stars are missing!!!

coords_data = []
line_num = 0
temp = [0,1]

Star_Data = []
No_Data = []
Sectors = []

for line in open('nonBlazRRC.dat'):
    parse = line.split("\t")
    temp[0] = parse[0]
    temp[1] = parse[1]
    coords_data.append(temp[:])
    
for r in coords_data:
    coords = SkyCoord(float(r[0]),float(r[1]), unit="deg")
    sector_table = Tesscut.get_sectors(coordinates=coords)
    print(sector_table)
    #if(len(sector_table)==0):
        #No_Data.append(coords)
    #else:
        #Sectors.append(sector_table)
    
#print(Sectors)
    #Star_Data.append(eleanor.multi_sectors(coords = coords, sectors = all))
    


#coords = SkyCoord(ra=68.959732, dec=-64.02704, unit=(u.deg, u.deg))

#star = eleanor.Source(coords=coords)
#data = eleanor.TargetData(star, height=15, width=15,uk )

#window =TK()
#window.title("Some Title")

