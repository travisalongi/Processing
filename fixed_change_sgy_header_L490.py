"""
Updated Aug. 2021
Updated Jan. 2023 - Fixed errors lines with gaps ~5 lines.

*** This is done on raw pre-stack segy that has header information and navigation 

Update segy header with navigation file information that contains time1 & lat/lon
processed all lines that are full lines, not turns *b.segy

@author: Travis Alongi
"""

import os, time, glob, utm, datetime
import pandas as pd
import numpy as np
import segyio as sgy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from obspy.core import UTCDateTime

# set some function definitions
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


# Import navigation (gps) data
nav = pd.read_csv(
    "l-4-90-sc.410",
    sep="\s+",
    skiprows=16,
    index_col=False,
    names=["datetime", "lat", "lon", "line", "?", "shot", "cdp"],
)

# Format datetime to a usable UTCDateTime format
times = []
for t in nav.datetime.astype(str).values:
    time1 = UTCDateTime(
        year=int(t[:4]),
        julday=int(t[4:7]),
        hour=int(t[7:9]),
        minute=int(t[9:11]),
        second=int(t[11:13]),
    )
    times.append(time1)
nav["utc_datetime"] = times

# Remove fraction of second from datetime
nav["datetime"] = [int(t[:-1]) for t in nav.datetime.astype(str).values]

# Create new columns of line numbers
nav_line = nav.line.to_list()
nav_num = [l.split("-")[-1] for l in nav_line]
nav["line_num"] = nav_num

# Get file names / paths
data_dir = "legacy_cmp_srt"
files = ["105", "107",  "113", "125"]
files = [data_dir + "/l4srt" + f + ".sgy" for f in files]

# Set segy byte information
sgt = sgy.TraceField
byte_yr = sgt.YearDataRecorded
byte_jday = sgt.DayOfYear
byte_hr = sgt.HourOfDay
byte_min = sgt.MinuteOfHour
byte_sec = sgt.SecondOfMinute

# Loop through and fix lines
plt.figure()
cmap = plt.cm.viridis(np.linspace(0, 1, len(files)))  # set color map for lines
for j, file in enumerate(files[:]):
    print("Processing file {}".format(file))
    print("Starting file at {}".format(str(datetime.datetime.now())))
    tic = time.time()
    line_number = file.split("/")[-1].split(".")[0].split("t")[-1]
    f = sgy.open(file, "r+", ignore_geometry=True)

    h = f.header  # All trace headers
    n = len(h)  # Total number of traces

    # Navigation subset
    nav_line = nav[nav.line_num == line_number]

    yr, jday, hr_list, mn_list, sec_list = [], [], [], [], []
    for i in range(n):
        y = str(h[i][byte_yr])
        jd = str(int((h[i][byte_jday])))
        hr = str(h[i][byte_hr])
        mn = str(h[i][byte_min])
        ss = str(h[i][byte_sec])

        # Format the header times
        if int(hr) < 10:
            hr = "0" + hr
        if int(mn) < 10:
            mn = "0" + mn
        if len(ss) > 3:
            ss = ss[:2]
        else:
            ss = "0" + ss[0]

        if int(hr) > 24:
            print("hour error. index ", i)
        if int(mn) > 60:
            print("min errror; index ", i)

        # Issue with julian day on line 107 {Hackk fix}
        if line_number == '107':
            yr.append(y)
            jday.append(131)
            hr_list.append(hr)
            mn_list.append(mn)
            sec_list.append(ss)

        else:
            yr.append(y)
            jday.append(jd)
            hr_list.append(hr)
            mn_list.append(mn)
            sec_list.append(ss)

    date_trace, date_trace_utc, err_mask = [], [], []
    for y, jd, hr1, m, s in zip(yr, jday, hr_list, mn_list, sec_list):
        # This is specifically for this survey: time1 not always correct
        if int(y) == 0:
            err_mask.append(False)
            # Place holders for traces w/o time1 in headers
            date_trace.append("0")
            date_trace_utc.append(UTCDateTime(1970, 1, 1))
        else:
            date_trace.append("19" + y + str(jd) + hr1 + m + s)
            date_trace_utc.append(
                UTCDateTime(
                    year=1990,
                    julday=int(jd),
                    hour=int(hr1),
                    minute=int(m),
                    second=int(s),
                )
            )
            err_mask.append(True)
    date_trace = np.array(date_trace)
    date_trace_utc = np.array(date_trace_utc)

    # Interpolate lat/lon as a function of trace time1
    f_lon, f_lat = parametric_interpolation(
        nav_line.lon, nav_line.lat, nav_line.datetime
    )

    dt_clean = np.where(date_trace != "0", date_trace, np.nan)
    trace_lon = f_lon(dt_clean)
    trace_lat = f_lat(dt_clean)

    # Populate src x&y headers
    for i, (lat, lon) in enumerate(zip(trace_lat, trace_lon)):
        if err_mask[i] == True:
            coords = utm.from_latlon(lat, lon)
            easting = remove_decimal([coords[0]])
            easting = int(easting[0])

            northing = remove_decimal([coords[1]])
            northing = int(northing[0])

            # Assign source locations
            h[i][sgy.TraceField.SourceX] = easting
            h[i][sgy.TraceField.SourceY] = northing
            h[i][
                sgy.TraceField.SourceGroupScalar
            ] = -100  # scaler that recomputes to account for moving decimal

        else:
            h[i][sgy.TraceField.SourceX] = 0
            h[i][sgy.TraceField.GroupX] = 0
            h[i][sgy.TraceField.INLINE_3D] = 0

    # Sanity checks
    plt.plot(nav_line.lon, nav_line.lat, "r.", label="Nav-file")
    plt.plot(
        trace_lon,
        trace_lat,
        ".",
        markeredgecolor=cmap[j],
        label="Header Value %s" % line_number,
        color=cmap[j],
    )
    plt.legend()

    toc = time.time()
    print(file, " run time1 = ", toc - tic, "\n")

    f.close
plt.show()
