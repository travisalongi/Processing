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
byte_srcx = sgt.SourceX
byte_srcy = sgt.SourceY

# Loop through and fix lines
for j, file in enumerate(files[7:]):
    print("Processing file {}".format(file))
    print("Starting file at {}".format(str(datetime.datetime.now())))
    tic = time.time()
    line_number = file.split("/")[-1].split(".")[0].split("t")[-1][:3]
    f = sgy.open(file, "r+", ignore_geometry=True)

    h = f.header    # All trace headers
    n = len(h)      # Total number of traces

    # Navigation subset
    nav_line = nav[nav.line_num == line_number]
    print(nav_line.shot.min(), nav_line.shot.max())
    
    ffids = np.zeros([n])
    for i in range(n):
        ffid = h[i][byte_ffid]
        ffids[i] = ffid
        if ffid in nav_line.shot.values:
    
            # Match to ffid in header to 410 file
            nav_trace = nav_line[nav_line.shot == ffid]
        
            # Normal condition navigation exists for ffid
            if nav_trace.size > 0:
                coords = utm.from_latlon(nav_trace.lat.values, nav_trace.lon.values)
                easting, northing = (
                    remove_decimal(coords[0])[0],
                    remove_decimal(coords[1])[0],
                )
        
                h[i][byte_srcx] = easting
                h[i][byte_srcy] = northing
        
                # Scaler that recomputes to account for moving decimal
                h[i][sgy.TraceField.SourceGroupScalar] = -100
        
            # Fix for when shot point (ffid is missing from 410)
            else:
                print("No Navigation for FFID %i" % ffid)
        
                # Take the mean of the shot loc. before and after the missing data
                nav_trace = nav_line[
                    (nav_line.shot <= ffid + 1) & (nav_line.shot >= ffid - 1)
                ]
                coords = utm.from_latlon(nav_trace.lat.mean(), nav_trace.lon.mean())
                easting, northing = (
                    remove_decimal([coords[0]])[0],
                    remove_decimal([coords[1]])[0],
                )
        
                h[i][byte_srcx] = easting
                h[i][byte_srcy] = northing
        
                # Scaler that recomputes to account for moving decimal
                h[i][sgy.TraceField.SourceGroupScalar] = -100
        
    toc = time.time()
    print(file, " run time1 = ", toc - tic, "\n")
    f.close
