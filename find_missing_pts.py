#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 10:17:19 2020

Imports navigation from csv, populates missing fields
Makes plots of each line with given points and extrapolates missing points
Saves location data to new csv

@author: travisalongi
"""

from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def parametric_interpolation(x ,y, t, type = 'linear'):
    """takes arrays of x, y, t parameterizes on t and returns interpolation
    of f(x(t)) and f(y(t)) as a tuple"""
    fx_t = interp1d(t,x, fill_value='extrapolate')
    fy_t = interp1d(t,y, fill_value='extrapolate')    
    return fx_t, fy_t

def make_plot(lat, lon, new_lat, new_lon, save_fig = False, fname = ''):
    final_point = [str(new_lon)[:6] + ' ' + str(new_lat)[:8]  ] #format string
    
    fig = plt.figure()
    buffer = 0.025 #for moving text string
    plt.plot(lat, lon, 'k:o', label = 'Given Pts.')
    plt.plot(new_lat, new_lon,  'ro', label = 'Extrapolated Pt.')
    plt.text( new_lat - buffer, new_lon + buffer, final_point, color = 'red')
    
    plt.legend(frameon=False, loc = 'upper left')
    
    plt.title(fname)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.axis('equal')
    plt.tight_layout()
    
    if save_fig == True:
        plt.savefig(fname + '.pdf')
    

    
#import data
df = pd.read_csv('B-1234-72-SC_line_nav.csv')   
df_new = df.copy()

missing_data = df.lat.isna() #find missing data points

df_missing = df[missing_data] #make new dataframes
df_filled = df[~missing_data]


for index, row in df_missing[:].iterrows():    
    line = row.line #number of line
    time = row.time
    df_line = df_filled[df_filled.line == line] #data frame of line
    
    print(line, len(df_line))
    
    if len(df_line) > 1: #needs 2 or more points to calculate
        x = df_line.lon
        y = df_line.lat
        t = df_line.time
        (fx_t, fy_t) = parametric_interpolation(x,y,t)
        
        x_new = fx_t(row.time)
        y_new = fy_t(row.time)
        
        df_new.at[index, 'lat'] = y_new
        df_new.at[index, 'lon'] = x_new   
        
        make_plot(x,y,x_new,y_new, save_fig=True, fname = line +'_'+ str(time))
        
    
#save file    
df_new.to_csv('B-1234-72-SC_line_nav_TA.csv', float_format = '%.5f', index = False)   