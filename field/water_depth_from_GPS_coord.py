#%%
import numpy as np
import pickle

import csv
from datetime import datetime
import pytz
import glob 
import scipy
import sys

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/')

from icewave.tools import weather

#%%

interp_H = weather.get_bathy_interpolator(disk='B:',year='2025')


# %% importer les données gps des acquisitions 
# de géophones et en déduire la bathymétrie à ces endroits
csv_file_coord = 'B:/General/Summary_geophone_lines/' + 'all_geophone_lines_coordinates.csv'

with open(csv_file_coord, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data = []
    for lines in csvFile:
        if count==0:
            header = lines
        else:
            data.append(lines)
        count+=1

print(data)    

list_H = []
for i in range(len(data)):
    H = float(weather.get_bathymetry_GPS((float(data[i][3]), float(data[i][2])), interp_H))
    list_H.append(H)
    # format : latitude, longitude
print(list_H)
# %%
# à continuer... (avoir aussi les heures des acquisitions pour trouver les marées...)
csv_file_times_UTC = 'B:/General/Summary_geophone_lines/' + 'all_geophone_lines_times_UTC.csv'

with open(csv_file_times_UTC, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data1 = []
    for lines in csvFile:
        if count==0:
            header1 = lines
        else:
            data1.append(lines)
        count+=1
print(header1)

print(data1) 

acq_numbers = []
for i in range(len(data1)):
    UTC_datetime_str = data1[i][2] # on prend l'instant où les géophones ont été allumés
    print(UTC_datetime_str)
    UTC_datetime = datetime.fromisoformat(UTC_datetime_str.replace("Z", "+00:00"))
    tide_height = weather.tide_from_datetime(UTC_datetime,disk = 'B:',year = '2025')
    print(tide_height)
# Parse ISO 8601 string
#dt = datetime.fromisoformat(dt_str.replace("Z", "+00:00"))
# %%
type(times_UTC_start[0])
# %%
