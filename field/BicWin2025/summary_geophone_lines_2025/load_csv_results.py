#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
from datetime import datetime

#%%
disk = 'B:'
csvfilepath = f"{disk}/General/Summary_geophone_lines/results_acquisitions.csv"


with open(csvfilepath, mode ='r')as file:
    csvFile = csv.reader(file, delimiter=',')
    count=0
    data_csv = []
    for lines in csvFile:
        if count==0:
            header = lines
        elif count>=1:
            data_csv.append(lines)
        count+=1

data_csv = np.array(data_csv)

print(header)

#%%
def convert_array1d_to_floatarray(arr):
    arr_new = []
    for i in range(len(arr)):
        if arr[i]=='':
            arr_new.append(np.nan)
        else:
            arr_new.append(float(arr[i]))
    return np.array(arr_new)


def convert_array1d_to_datetimearray(arr):
    arr_new = []
    for i in range(len(arr)):
        dt = datetime.strptime(arr[i], "%Y-%m-%dT%H:%M:%S.%fZ")
        arr_new.append(dt)
    return np.array(arr_new)   

# %%
%matplotlib qt

plt.figure()
plt.plot(convert_array1d_to_datetimearray(data_csv[:,5]), convert_array1d_to_floatarray(data_csv[:,7]), '.b')
plt.ylim(0, 7e9)
plt.show()
# %%
