#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 13:50:44 2020
Update segy header with navigation file information that contains time & lat/lon
processed all lines that are full lines, not turns *b.segy
@author: Travis Alongi
"""
import os 
import time
import glob
import pandas as pd
import numpy as np
import segyio as sgy
import conversion as utm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from obspy.core import UTCDateTime
from collections import Counter

#set some function definitions
def remove_decimal(number, places_after_decimal_to_retain = 2):
    """removes decimal place and truncates the new number"""
    n_list_strings = str(number).split('.')
    truncated_decimal = n_list_strings[1][:places_after_decimal_to_retain]
    new_number = int(n_list_strings[0] + truncated_decimal)
    return new_number

def parametric_interpolation(x ,y, t, type = 'linear'):
    """takes arrays of x, y, t parameterizes on t and returns interpolation
    of f(x(t)) and f(y(t)) as a tuple"""
    fx_t = interp1d(t,x, fill_value='extrapolate')
    fy_t = interp1d(t,y, fill_value='extrapolate')    
    return fx_t, fy_t
 
    
## import data
nav = pd.read_csv('/run/media/admin/mojo128/a100sc/a-1-00-sc.065_410_sparker_multichannel', 
                       sep = '\s+',
                       skiprows = 44,
                       index_col = False,
                       names = ['datetime', 'lat', 'lon', 'ffid', 'shot', 'cdp'])

# convert nav to utc time
nav['utc_datetime'] = [UTCDateTime(year = int(t[:4]), 
                          julday = int(t[4:7]), 
                          hour = int(t[7:9]), 
                          minute = int(t[9:11]), 
                          second = int(t[11:13])) for t in nav.datetime.astype(str).values]

 #remove fraction of second from datetime
nav['datetime'] = [int(t[:-1]) for t in nav.datetime.astype(str).values]


# switch directories get file names
data_dir = '/run/media/admin/mojo128/a100sc/a100sc_shots'
os.chdir(data_dir)
#files = glob.glob('*.sgy')
files = ['401.sgy', '830.sgy', '847_1.sgy', 'a1000248.sgy']

cmap = plt.cm.viridis(np.linspace(0,1,len(files))) #set color map for lines


#initiate plot (part of sanity check)
plt.figure(1,figsize = (12,12))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Line map for A100SC')

#%%

## loop through files, change header
for j, file in enumerate(files[-1:]):
    tic = time.time()
    line_number = file.split('.')[0][-4:]
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
#    yr = [str(h[i][byte_yr]) for i in range(n)]  
#    jday = [str(h[i][byte_jday]) for i in range(n)]
#    if np.sum(np.diff(np.array(jday).astype(int)) != 0) > 1:
#        print('prossible header error check julian dates')
#        break
    
    yr = []
    jday = [h[0][byte_jday]]
    hr_list = []
    mn_list = []
    sec_list = []
    rm_err = np.ones(n)
    for i in range(n):
        y = str(h[i][byte_yr])
        jd = str(int((h[i][byte_jday])))
        hr = str(h[i][byte_hr])
        mn = str(h[i][byte_min])
        ss = str(h[i][byte_sec])
        
        
        if i > 0:
            if np.abs(int(jd) - int(jday[0])) > 1:
                print('Skipping trace # ' + str(i) + ' in '  + file + ' a portion due to garbage header')
                rm_err[i] = 0
                continue
            
                      
        if int(hr) < 10:
            hr = '0' + hr
        if int(mn) < 10:
            mn = '0' + mn
        if int(ss) < 10:
            ss = '0' + ss
            
        if int(hr) > 24:
            print('hour error. index ', i)
        if int(mn) > 60:
            print('min errror; index ', i)

        yr.append(y)    
        jday.append(jd)
        hr_list.append(hr)
        mn_list.append(mn)
        sec_list.append(ss)


    date_trace_utc = []
    date_trace = []
    for y,jd,hr1,m,s in zip(yr,jday,hr_list,mn_list, sec_list):
        date_trace.append('200' + y +str(jd) + hr1 + m + s)
        date_trace_utc.append(UTCDateTime(year = 2000, 
                                  julday = int(jd), 
                                  hour = int(hr1), 
                                  minute = int(m), 
                                  second = int(s)))
      
#    date_trace = [int('200' + y + j + hr1 + m + s) for y,j,hr1,m,s in zip(yr,jday,hr_list,mn_list,sec_list)]    
#    date_trace_utc = [UTCDateTime(year = int('200' + y), 
#                                  julday = int(j), 
#                                  hour = int(hr2), 
#                                  minute = int(m), 
#                                  second = int(s)) for y,j,hr2,m,s in zip(yr,jday,hr_list,mn_list,sec_list)]

    
    ##trim nav file for trace times
    line_time_begin = date_trace_utc[0]
    line_time_end = date_trace_utc[-1]
    mask = (nav.utc_datetime >= line_time_begin) & (nav.utc_datetime <= line_time_end)
    
    #nav for this line
    nav_line = nav[mask]
    
    #interpolate lat/lon as a function of trace time
    f_lon, f_lat = parametric_interpolation(nav_line.lon, nav_line.lat, nav_line.datetime)
    trace_lon = f_lon(date_trace)
    trace_lat = f_lat(date_trace)
    
    #convert trace location to utm
    utm_coords = utm.from_latlon(trace_lat, trace_lon)
    
    #populate src x&y headers
    for i, (x,y) in enumerate(zip(utm_coords[0], utm_coords[1])):
        if i == rm_err[i]:
            easting = remove_decimal(x)
            northing = remove_decimal(y)
            
            h[i][sgy.TraceField.SourceX] = easting #assign source locations
            h[i][sgy.TraceField.SourceY] = northing
            h[i][sgy.TraceField.SourceGroupScalar] = -100 #scaler that recomputes to account for moving decimal
        
    #sanity checks
    plt.plot(nav_line.lon, nav_line.lat, 'k.')
    plt.plot(trace_lon, trace_lat, '.', color = cmap[j])
    plt.text(trace_lon[0], trace_lat[0], line_number)
    sample_rate_counts = Counter(np.diff(np.unique(date_trace_utc)))
    print('Line #',line_number,'**sample rates =' , sample_rate_counts)
    
    toc = time.time()
    print(file, ' run time = ',toc-tic, '\n')
    
    f.close
    
    

