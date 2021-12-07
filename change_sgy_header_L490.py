#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated Aug. 2021

*** This is done on raw pre-stack segy that has header information and navigation 

Update segy header with navigation file information that contains time & lat/lon
processed all lines that are full lines, not turns *b.segy

@author: Travis Alongi
"""

import os, time, glob, utm, datetime
import pandas as pd
import numpy as np
import segyio as sgy
import matplotlib.pyplot as plt
from scipy.stats import mode
from scipy.interpolate import interp1d
from obspy.core import UTCDateTime
from collections import Counter

#set some function definitions
def remove_decimal(numbers, places_after_decimal_to_retain = 2):
    """removes decimal place and truncates the new number"""
    new_numbers = []
    for number in numbers:
        n_list_strings = str(number).split('.')
        truncated_decimal = n_list_strings[1][:places_after_decimal_to_retain]
        new_number = int(n_list_strings[0] + truncated_decimal)
        new_numbers.append(new_number)
    return new_numbers

def parametric_interpolation(x ,y, t, type = 'linear'):
    """takes arrays of x, y, t parameterizes on t and returns interpolation
    of f(x(t)) and f(y(t)) as a tuple"""
    fx_t = interp1d(t,x, fill_value='extrapolate')
    fy_t = interp1d(t,y, fill_value='extrapolate')    
    return fx_t, fy_t


    
# Import data
nav = pd.read_csv('l-4-90-sc.410', 
                       sep = '\s+',
                       skiprows = 16,
                       index_col = False,
                       names = ['datetime', 'lat', 'lon', 'ffid', 'shot', 'cdp'])

# Format datetime to a usable UTCDateTime format
times = []
for t in nav.datetime.astype(str).values:
    time1 = UTCDateTime(year = int(t[:4]), 
                          julday = int(t[4:7]), 
                          hour = int(t[7:9]), 
                          minute = int(t[9:11]), 
                          second = int(t[11:13]))
    times.append(time1)
nav['utc_datetime'] = times


# Remove fraction of second from datetime
nav['datetime'] = [int(t[:-1]) for t in nav.datetime.astype(str).values]


# Switch directories get file names
data_dir = 'legacy_cmp_srt'
os.chdir(data_dir)
files = glob.glob('*.sgy')

cmap = plt.cm.viridis(np.linspace(0,1,len(files))) #set color map for lines

#%%

# Byte information
byte_yr = sgy.TraceField.YearDataRecorded
byte_jday = sgy.TraceField.DayOfYear
byte_hr = sgy.TraceField.HourOfDay
byte_min = sgy.TraceField.MinuteOfHour
byte_sec = sgy.TraceField.SecondOfMinute

## Loop through files changing headers
for j, file in enumerate(files):
# for j, file in enumerate(['l4srt101.sgy']):
    print('Processing file {}'.format(file))
    print('Starting file at {}'.format(str(datetime.datetime.now())))
    tic = time.time()
    line_number = file.split('.')[0].split('t')[1]
    f = sgy.open(file, 'r+', ignore_geometry = True)
    
    h = f.header # All trace headers
    n = len(h) # Total number of traces
    
    
    yr = []
    jday = []
    hr_list = []
    mn_list = []
    sec_list = []

    for i in range(n):

        y = str(h[i][byte_yr])
        jd = str(int((h[i][byte_jday])))
        hr = str(h[i][byte_hr])
        mn = str(h[i][byte_min])
        ss = str(h[i][byte_sec])
         
        # Format the header times              
        if int(hr) < 10:
            hr = '0' + hr
        if int(mn) < 10:
            mn = '0' + mn
        # if int(ss) < 10:
        #     ss = '0' + ss
        if len(ss) > 3:
            ss = ss[:2]
        else:
            ss = '0' + ss[0]
            
        if int(hr) > 24:
            print('hour error. index ', i)
        if int(mn) > 60:
            print('min errror; index ', i)

        yr.append(y)    
        jday.append(jd)
        hr_list.append(hr)
        mn_list.append(mn)
        sec_list.append(ss)
     
    ffid_list, bin_x, bin_y, cdp = [], [], [], []    
    for i in range(n):
        ffid_list.append(h[i][9])
        bin_x.append(h[i][181])
        bin_y.append(h[i][185])
        cdp.append(h[i][sgy.TraceField.CDP])
        
    
    jd_mode = int(mode(jday)[0])
    date_trace_utc = []
    date_trace = []
    err_mask = [] # Boolean where bad traces are False
    for y,jd,hr1,m,s in zip(yr,jday,hr_list,mn_list, sec_list):
        
        # This is specifically for this survey because time info not always correct in header
        if int(y) == 0 or (int(jd) <= jd_mode - 1 or int(jd) >= jd_mode + 1): 
            # Boolean
            err_mask.append(False)
            # Place holders for traces w/o time in headers
            date_trace.append('0')
            date_trace_utc.append(UTCDateTime(1970,1,1))
        else:
            date_trace.append('19' + y +str(jd) + hr1 + m + s)
            date_trace_utc.append(UTCDateTime(year = 1990, 
                                      julday = int(jd), 
                                      hour = int(hr1), 
                                      minute = int(m), 
                                      second = int(s)))
            err_mask.append(True)
    date_trace = np.array(date_trace)
    date_trace_utc = np.array(date_trace_utc)
      
    
    # Trim nav file between the time of segy file / line
    line_time_begin = (date_trace_utc)[err_mask][0]
    line_time_end = (date_trace_utc)[err_mask][-1]
    mask = (nav.utc_datetime >= line_time_begin) & (nav.utc_datetime <= line_time_end)
    
    # Nav for this line
    nav_line = nav[mask]
    
    # Interpolate lat/lon as a function of trace time
    f_lon, f_lat = parametric_interpolation(nav_line.lon, nav_line.lat, nav_line.datetime)
    trace_lon = f_lon(date_trace[err_mask])
    trace_lat = f_lat(date_trace[err_mask])
  
    # Populate src x&y headers
    # for i, (x,y) in enumerate(zip(utm_coords[0], utm_coords[1])):
    for i, (lat,lon) in enumerate(zip(trace_lat, trace_lon)):
        # print(i,lat,lon)
        if err_mask[i] == True:
            
            coords = utm.from_latlon(lat, lon)
            
            easting = remove_decimal(coords[0])
            northing = remove_decimal(coords[1])
            
            h[i][sgy.TraceField.SourceX] = easting #assign source locations
            h[i][sgy.TraceField.SourceY] = northing
            h[i][sgy.TraceField.SourceGroupScalar] = -100 #scaler that recomputes to account for moving decimal
        
        else:
            h[i][sgy.TraceField.SourceX] = 0
            h[i][sgy.TraceField.GroupX] = 0
            h[i][sgy.TraceField.INLINE_3D] = 0
            
#    # Sanity checks
#    plt.figure()
#    plt.plot(nav_line.lon, nav_line.lat, 'k.')
#    plt.plot(trace_lon, trace_lat, '.', color = cmap[j])
#    plt.text(trace_lon[0], trace_lat[0], line_number)
#    print('Line #',line_number)
    
    toc = time.time()
    print(file, ' run time = ',toc-tic, '\n')
    
    f.close
