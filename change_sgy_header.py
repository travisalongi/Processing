#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:50:44 2020
Update segy header with navigation file information that contains time & lat/lon
@author: talongi
"""
import os 
import time
import glob
import pandas as pd
import numpy as np
import segyio as sgy
import utm
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime

def remove_decimal(number, places_after_decimal_to_retain = 2):
    """removes decimal place and truncates the new number"""
    n_list_strings = str(number).split('.')
    truncated_decimal = n_list_strings[1][:places_after_decimal_to_retain]
    new_number = int(n_list_strings[0] + truncated_decimal)
    return new_number
    
## import data
nav = pd.read_csv('E:\\O199SC\\o199sc\\nav\\o-1-99-sc.061', 
                       sep = ' ',
                       skiprows = 10,
                       names = ['datetime', 'lat', 'lon'])

data_dir = 'E:\\O199SC\\airgun'
os.chdir(data_dir)
files = glob.glob('*.sgy')
cmap = plt.cm.viridis(np.linspace(0,1,len(files)))


## loop through files, change header
for j, file in enumerate(files):
    f = sgy.open(file, 'r+', ignore_geometry = True)
    
    h = f.header #all trace headers
    n = len(h) #total number of traces
    
    ## handle date formating
    # byte information
    byte_yr = sgy.TraceField.YearDataRecorded
    byte_jday = sgy.TraceField.DayOfYear
    byte_hr = sgy.TraceField.HourOfDay
    byte_min = sgy.TraceField.MinuteOfHour
    byte_sec = sgy.TraceField.SecondOfMinute
    
    # date information from trace headers
    yr = [str(h[i][byte_yr]) for i in range(n)]  
    jday = [str(h[i][byte_jday]) for i in range(n)]
    
    hr_list = []
    mn_list = []
    sec_list = []
    for i in range(n):
        hr = str(h[i][byte_hr])
        mn = str(h[i][byte_min])
        ss = str(h[i][byte_sec])
        
        if int(hr) < 10:
            hr = '0' + hr
        if int(mn) < 10:
            mn = '0' + mn
        if int(ss) < 10:
            ss = '0' + ss
        
        hr_list.append(hr)
        mn_list.append(mn)
        sec_list.append(ss)
                
    date_trace = [int('19' + y + j + h + m + s) for y,j,h,m,s in zip(yr,jday,hr_list,mn_list,sec_list)]    
    date_nav = [int(str(d)[:-1]) for d in nav.datetime]
    
    
    ## compare date/time from gps navigation to traces and populate source x&y location 
    tic = time.time()
    for i, date in enumerate(date_trace):
        index_closest_time = np.abs(date - np.array(date_nav)).argmin() #gets index of closest time between trace and navigation
        lat_near = nav.lat[index_closest_time]
        lon_near = nav.lon[index_closest_time]
        utm_coords = utm.from_latlon(lat_near, lon_near) #calculate utm coord b/c segy headers do not take floats
        
        easting = remove_decimal(utm_coords[0]) #remove decimal using function above
        northing = remove_decimal(utm_coords[1])
        
        h[i][sgy.TraceField.SourceX] = easting #assign source locations
        h[i][sgy.TraceField.SourceY] = northing
        h[i][sgy.TraceField.SourceGroupScalar] = -100 #scaler that recomputes to account for moving decimal
    toc = time.time()
    print(file, ' run time = ',toc-tic)
    
    byte_number = sgy.TraceField.SourceX
    src_x = [h[i][byte_number] for i in range(n)]
    src_y = [h[i][77] for i in range(n)]
    
    plt.figure(1)
    plt.plot(src_x, src_y, label = file)
    plt.axis('equal')
    plt.legend()
    
    f.close
