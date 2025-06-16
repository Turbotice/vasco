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
import sys
#%%
sys.path.append('C:/Users/Vasco/Zanchi/Documents/git_turbotice/')
from vasco.tools import open_csv
# %% creation d'un dictionnaire pour ranger les résultats obtenus

# inputs :
(dict_lists , list_categories , list_types) = open_csv.open_csv_table_experiments(file_path='R:/Gre25/Summary/dispersion_relations/inputs_RDD.csv')




print(list_categories)
print(list_types)


general_folder = f'R:/Gre25/Data/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'


dict_results = {}
for i in range(len(dict_lists['\ufefflist_dates'])):
    dict_results[dict_lists['\ufefflist_dates'][i]] = {}
for i in range(len(dict_lists['\ufefflist_dates'])):
    dict_results[dict_lists['\ufefflist_dates'][i]]['acq_'+str(dict_lists['list_acq_num'][i])] = {}    
for i in range(len(dict_lists['\ufefflist_dates'])):
    dict_results[dict_lists['\ufefflist_dates'][i]]['acq_'+str(dict_lists['list_acq_num'][i])]['camera_'+str(dict_lists['list_camera_SN'][i])] = {}

for i in range(len(dict_lists['\ufefflist_dates'])):
    dd = dict_results[dict_lists['\ufefflist_dates'][i]]['acq_'+str(dict_lists['list_acq_num'][i])]['camera_'+str(dict_lists['list_camera_SN'][i])]
    dd['W'] = dict_lists['list_W'][i]
    dd['Dt'] = dict_lists['list_Dt'][i]
    #dd['tab_v_phase'] , dd['tab_v_phase_err'] = routine_vitesse_phase(f_exc=tab_f_exc[j],freq_acq=tab_freq_acq[j],general_folder=general_folder,W=dd['W'],Dt=dd['Dt'],index_profile_line=dict_lists['list_indices_profile_line'][i],xlim_fit=(dict_lists['list_xinf'][i],dict_lists['list_xsup'][i]),dcm=dict_lists['list_dcm'][i],dpx=dict_lists['list_dpx'])
    
    general_folder = f'R:/Gre25/Data/{dict_lists['\ufefflist_dates'][i][-4:]}/cameras/manip_relation_dispersion/Acquisition_{str(dict_lists['list_acq_num'][i])}/camera_{dict_lists['list_camera_SN'][i]}/'
    print(general_folder)
    try:
        tab_f_exc,tab_freq_acq = list_all_freq(general_folder=general_folder)
        dd['tab_f_exc'] = tab_f_exc
        dd['tab_freq_acq'] = tab_freq_acq
        tab_v_phase = np.zeros(len(tab_f_exc))
        tab_v_phase_err = np.zeros(len(tab_f_exc)) 
        for j in range(len(tab_f_exc)):
            tab_v_phase[j],tab_v_phase_err[j] = routine_vitesse_phase(f_exc=tab_f_exc[j],freq_acq=tab_freq_acq[j],general_folder=general_folder,W=dd['W'],Dt=dd['Dt'],index_profile_line=dict_lists['list_indices_profile_line'][i],xlim_fit=(dict_lists['x_inf_fit'][i],dict_lists['x_sup_fit'][i]),dcm=dict_lists['list_dcm'][i],dpx=dict_lists['list_dpx'][i],plot=False)
        dd['tab_v_phase'] = tab_v_phase
        dd['tab_v_phase_err'] = tab_v_phase_err
    except:
        pass

    # enlever le try except (il faut avoir tous les fichiers piv prêts...)


#%%

for date in dict_results:
    for acq in dict_results[date]:
        for cam in dict_results[date][acq]:
            ddd = dict_results[date][acq][cam]
#            print(ddd)
            try:
                plt.plot(ddd['tab_f_exc'],ddd['tab_v_phase'],'o')
                plt.ylim(0,50)
            except:
                pass




























################################################################
# FIN DE LA PARTIE OÙ ON RUN POUR PLEIN D'ACQUISITION ##########
################################################################
#%%
plt.figure(figsize=(15,10))
for idx_routine in range(len(list_dates)):
    tab_f_exc = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_f_exc']
    tab_freq_acq = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_freq_acq']
    tab_v_phase = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase']
    tab_v_phase_err = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase_err']
    #print(tab_v_phase)
    mask = (tab_v_phase_err < 1/20*tab_v_phase) & (tab_v_phase>0)
    plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='o',linestyle='',label=list_dates[idx_routine])
    #plt.errorbar(tab_f_exc,tab_v_phase,tab_v_phase_err,marker='.',linestyle='')
plt.ylim(1,40)
plt.xlabel('frequency (Hz)')
plt.ylabel('phase velocity (m/s)')
plt.legend()
#plt.loglog()
plt.show()
#%% fitter en regime flexion pour avoir modul de flexion
def compute_omega(k,D):
    return np.sqrt(((D/rho)*k**5 + (T/rho)*k**3 + g*k)*np.tanh(k*H))
def compute_freq(k,D):
    return (1/(2*np.pi))*compute_omega(k,D)





"""
e = 2e-3 
H = 14.5e-2 - e 
E = 9e9
nu = 0.4
D_value = (E*(e**3)/(12*(1-nu**2)))
rho = 1e3
T = 0
lambdamin = 5e-2
lambdamax = 1
g = 9.81
tab_lambda = np.linspace(lambdamin,lambdamax,10000)
"""
########################## ICI CHOISIR "L'INDICE ROUTINE" A AFFICHER ET FITTER ##############
list_idx_routine=[0,1,15,16]
#############################################################################################

plt.figure(figsize=(15,10))

for idx_routine in list_idx_routine:
    tab_f_exc = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_f_exc']
    tab_freq_acq = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_freq_acq']
    tab_v_phase = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase']
    tab_v_phase_err = dico[list_dates[idx_routine]]['acq'+str(list_acq_num[idx_routine])]['camera_'+list_camera_SN[idx_routine]][f'W{str(list_W[idx_routine])}_Dt{str(list_Dt[idx_routine])}'][f'index_profile_line{str(list_indices_profile_line[idx_routine])}'][f'_xlim_fit{str(list_xlim_fit[idx_routine])}']['tab_v_phase_err']


    mask = (tab_v_phase_err < 1/10*tab_v_phase) & (tab_v_phase>0)
    plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='o',linestyle='',label=list_dates[idx_routine]+' ; acq'+str(list_acq_num[idx_routine])+' ; cam '+str(list_camera_SN[idx_routine]))


    popt,pcov = curve_fit(v_phase_flexural,tab_f_exc[mask],tab_v_phase[mask],sigma=tab_v_phase_err[mask])
    #plt.figure()
    tab_f_fit = np.linspace(10,200,1000)

    plt.plot(tab_f_fit,v_phase_flexural(tab_f_fit,popt[0]),label='Fit (if only flexural) : D = ('+str(np.round(popt[0],decimals=2))+ '+-'+ str(np.round(np.sqrt(pcov[0][0]),decimals=2)) + ') J')

    #plt.errorbar(np.array(dico_dataset['tab_f_exc'])[mask],np.array(dico_dataset['v_phase'])[mask],np.array(dico_dataset['v_phase_err'])[mask],marker='o')
    #plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='.',linestyle='')
    plt.ylim(0,50)
    plt.xlim(0,300)
    plt.xlabel('frequency (Hz)',fontsize=15)
    plt.ylabel('$v_{\phi}$',fontsize=15)
plt.legend()
plt.show()

