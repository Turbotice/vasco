#%% import libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
#import fitutils as fu
from scipy.signal import find_peaks
import os
import pickle
import csv
import re
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.interpolate import LinearNDInterpolator
#matplotlib.use('TkAgg')
# %%
#%matplotlib widget
%matplotlib qt
#%%


# %% load data
W = 64
#Dt = 1 # pour ce cas pas de Dt car on compare tout par rapport à une même image de reference
i0 = 7100
N = 9400
refimg = 0

date = '0522'
name_frac_file = 'recherche_freq_puis_frac'
#camera_SN = '22458101'

#f_exc = 0.94
freq_acq = 20

frame_frac = 9276
#computer = 'adour'

system_loc = 'windows_server'

#if computer=='DellVasco':
#    general_folder = f'K:/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#elif computer=='Leyre':
#    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_fracture/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

if system_loc=='linux_server':
    general_folder = f'/media/turbots/GreDisk/Gre25/Data/{date}/cameras/frac/image_sequence/'
elif system_loc=='windows_server':
    general_folder = f'R:/Gre25/Data/PIV_results/{date}_frac_{name_frac_file}/'


path2data = general_folder

matfile = f'{path2data}PIV_processed_3passages_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'

#matfile = f'{path2data}PIV_processed_i0{i0}_N{N}_W{W}_refimg{refimg}.mat'


from scipy.io import loadmat

mat_dict = loadmat(matfile)


u_original = mat_dict['u_original'][:,0]
v_original = mat_dict['v_original'][:,0]

def reshape_array(arr):
    array_new = np.zeros((len(arr),arr[0].shape[0],arr[0].shape[1]))
    for i in range(len(arr)):
        array_new[i,:,:] = arr[i]
    return array_new

u = reshape_array(u_original)
v = reshape_array(v_original)

u = -u
v = -v

#%%

yind = 4 # il vaut mieux utiliser les coordonnées en pixels

plt.figure()
plt.imshow(v[1500],vmin=-5,vmax=5)
plt.colorbar()
plt.show()


plt.figure()
for i in range(2170,2180): # fenetre temporelle où il y a la fracture
    plt.plot(v[i,yind,:],label=str(i+i0))
plt.legend()
plt.show()


# %% courbure de la plaque juste avant la fracture ?

xvals = np.arange(v.shape[2]) * 0.075 * W/2 * 1e-2
v_converted_meters = v * 0.075 * 1e-2

ind_inf_fit = 24
ind_sup_fit = 35
print(xvals)
fit_params = np.polyfit(xvals[ind_inf_fit:ind_sup_fit],v_converted_meters[frame_frac-i0,yind,ind_inf_fit:ind_sup_fit],2)
print(fit_params)

plt.figure()
plt.plot(xvals,v_converted_meters[frame_frac-i0,yind,:],label='frame '+str(i+i0))
plt.plot(xvals[ind_inf_fit:ind_sup_fit],(fit_params[0]*xvals**2 + fit_params[1]*xvals + fit_params[2])[ind_inf_fit:ind_sup_fit],label='fit : $\kappa$ = '+str(np.round(2*fit_params[0],3)))
plt.legend()
plt.show()

kappa_c = 2*fit_params[0]
# %%
