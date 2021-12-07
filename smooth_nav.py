#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:14:28 2021

This code takes the group XY locations of each trace in a line.
The x/y trace locations are then smoothed with a running mean over a window
The smoothed x/y are written to same segy file in CDP X/Y header locations

@author: talongi
"""

import os 
import time
import glob
import numpy as np
import pandas as pd
import segyio as sgy
import matplotlib.pyplot as plt
import utm

window = 10 #moving number of traces


# switch directories get file names
data_dir = '/home/talongi/Gypsy/Processing/Processing/O199SC_processed/3ms_gap'
save_dir = '/home/talongi/Gypsy/Processing/Processing/O199SC_processed/3ms_gap_smooth_nav/'
os.chdir(data_dir)
files = glob.glob('*.segy')

def remove_decimal(number, places_after_decimal_to_retain = 2):
    """removes decimal place and truncates the new number"""
    n_list_strings = str(number).split('.')
    truncated_decimal = n_list_strings[1][:places_after_decimal_to_retain]
    new_number = int(n_list_strings[0] + truncated_decimal)
    return new_number


## loop through files, change header
for j, file in enumerate(files[:]):
    tic = time.time()
    line_number = file.split('.')[0][-4:]
    f = sgy.open(file, 'r+', ignore_geometry = True)
    
    h = f.header #all trace headers
    n = len(h) #total number of traces

    # get group x & y for each trace
    xs, ys = [], []
    for i in h: # loop through all the traces
        x = i[sgy.TraceField.SourceX]
        y = i[sgy.TraceField.SourceY]
        xs.append(x)
        ys.append(y)
        
    plt.figure()
    plt.plot(xs,ys,'-o')
    
    # put into pd series to smooth
    X = pd.Series(xs)
    Y = pd.Series(ys)
    
    # tried rolling median and rolling mean -- mean appears to be more accurate to start and stop of line
    Xr = X.rolling(window, min_periods = 1, center = True, win_type = 'hamming').mean()
    Yr = Y.rolling(window, min_periods = 1, center = True, win_type = 'hamming').mean()
        
    # repopulate headers with smoothed navigation
    for j, (xn,yn) in enumerate(zip(Xr, Yr)):
        easting = remove_decimal(xn, places_after_decimal_to_retain=0)
        northing = remove_decimal(yn, places_after_decimal_to_retain=0)
        
        h[j][sgy.TraceField.CDP_X] = easting
        h[j][sgy.TraceField.CDP_Y] = northing
        
    toc = time.time()
    print('%s processed in %i seconds' % (file, toc - tic))
    f.close()
        
    