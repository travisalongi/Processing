""" Uses shots from 410 file to populate segy headers using utm 
shotpoint locations.

Updated Aug. 2021
Updated Jan. 2023 - Fixed errors lines with gaps ~5 lines.
Updated Mar. 2023 - Rewrite to match on ffids instead of interpolating the time

This is done on raw pre-stack segy that has header information and navigation 

@author: Travis Alongi
"""

import os, time, glob, utm, datetime
import pandas as pd
import numpy as np
import segyio as sgy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Set some function definitions
def remove_decimal(numbers, places_after_decimal_to_retain=2):
    """Removes decimal place and truncates the new number."""
    new_numbers = []
    for number in numbers:
        n_list_strings = str(number).split(".")
        truncated_decimal = n_list_strings[1][:places_after_decimal_to_retain]
        new_number = int(n_list_strings[0] + truncated_decimal)
        new_numbers.append(new_number)
    return new_numbers

def parametric_interpolation(x, y, t, type="linear"):
    """Takes arrays of x, y, t parameterizes on t and returns interpolation
    of f(x(t)) and f(y(t)) as a tuple"""
    fx_t = interp1d(t, x, fill_value="extrapolate")
    fy_t = interp1d(t, y, fill_value="extrapolate")
    return fx_t, fy_t

# Import navigation (GPS) data
nav = pd.read_csv(
    "l-4-90-sc.410",
    sep="\s+",
    skiprows=16,
    index_col=False,
    names=["datetime", "lat", "lon", "line", "?", "shot", "cdp"],
)


# Add line numbers to nav data frame
nav_line_list = nav.line.to_list()
nav_num = [l.split("-")[-1] for l in nav_line_list]
nav["line_num"] = nav_num

# Get file names / paths
data_dir = "legacy_cmp_srt"
files = glob.glob(data_dir + "/*.sgy")

# Set segy byte information
sgt = sgy.TraceField
byte_ffid = sgt.FieldRecord
byte_cdp = sgt.CDP
byte_srcx = sgt.SourceX
byte_srcy = sgt.SourceY

# Loop through and fix lines
plt.figure()
for j, file in enumerate(files[2:27]):
    print("Processing file {}".format(file))
    print("Starting file at {}".format(str(datetime.datetime.now())))
    tic = time.time()
    line_number = file.split("/")[-1].split(".")[0].split("t")[-1][:3]


    with sgy.open(file, "r+", ignore_geometry=True) as f:

        h = f.header    # All trace headers
        n = len(h)      # Total number of traces

        # Navigation subset
        nav_line = nav[nav.line_num == line_number]
        print(nav_line.cdp.min(), nav_line.cdp.max(), n, str(datetime.datetime.now()))

        # Parameterize navigation - to account for odd cpds in the segy
        fx, fy = parametric_interpolation(nav_line.lon, nav_line.lat, nav_line.cdp)
        
        ffids = np.zeros([n])
        cdps = np.zeros([n])
        xs, ys = np.zeros([n]), np.zeros([n])
        for i in range(n):

            ffid = h[i][byte_ffid]
            ffids[i] = ffid
            cdp = h[i][byte_cdp]
            cdps[i] = cdp

            # Interpolate the lat/lon from the cdp info
            mod_lon, mod_lat = fx(cdp), fy(cdp)
      
            # Normal condition navigation exists for ffid
            coords = utm.from_latlon(mod_lat, mod_lon)
            easting, northing = (
                remove_decimal([coords[0]])[0],
                remove_decimal([coords[1]])[0],
            )
  
            h[i][byte_srcx] = easting
            h[i][byte_srcy] = northing
            xs[i] = easting
            ys[i] = northing
  
            # Scaler that recomputes to account for moving decimal
            h[i][sgy.TraceField.SourceGroupScalar] = -100

                #          
                # # Fix for when shot point (ffid is missing from 410)
                # else:
                #     print("No Navigation for CDP %i" % cdp)
                #          
                #     # Take the mean of the shot loc. before and after the missing data
                #     nav_trace = nav_line[
                #         (nav_line.shot <= ffid + 1) & (nav_line.shot >= ffid - 1)
                #     ]
                #     coords = utm.from_latlon(nav_trace.lat.mean(), nav_trace.lon.mean())
                #     easting, northing = (
                #         remove_decimal([coords[0]])[0],
                #         remove_decimal([coords[1]])[0],
                #     )
                #          
                #     h[i][byte_srcx] = easting
                #     h[i][byte_srcy] = northing
                #          
                #     # Scaler that recomputes to account for moving decimal
                #     h[i][sgy.TraceField.SourceGroupScalar] = -100
          
    toc = time.time()
    print(file, " run time1 = ", toc - tic, "\n")
