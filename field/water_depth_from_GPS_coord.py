#%%
import numpy as np
import pickle

import csv
from datetime import datetime
import pytz
import glob 
import scipy
import sys

sys.path.append('C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave')

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
