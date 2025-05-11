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
#matplotlib.use('TkAgg')
# %% def fonctions
def savedict(filename,dico):
    with open(filename+'.pickle', 'wb') as handle:
        pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)
def loaddict(filename):
    with open(filename+'.pickle', 'rb') as handle:
        dico = pickle.load(handle)
    return dico
def extract_data_mat_file(path_mat_file):
    with h5py.File(path_mat_file, 'r') as file:
        data = file['data'][:]  # Accessing 'data' variable
    return data

def load_complex_field(path_mat_file):
    matrix = extract_data_mat_file(path_mat_file)
    # Initialize an empty matrix to store complex numbers
    complex_matrix = np.empty(matrix.shape, dtype=complex)
    
    # Iterate through each element in the matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            # Get the tuple from the original matrix
            a, b = matrix[i, j]
            
            # Create a complex number from the tuple and assign it to the corresponding position
            complex_matrix[i, j] = complex(a, b)
    return complex_matrix    

def load_csv_input_RDD(file_path = '/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/inputs_RDD.csv',skipfirst2rows = True):
    # Dictionaries to store categories and lists
    list_categories = []
    dict_lists = {}
    list_types = []

    # Open and read the CSV file
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        count = 0

        for row in reader:
            if count == 0 and skipfirst2rows:  # Handle the header row
                list_categories = [row[i] for i in range(len(row))]  # Map index to column name
                for i in range(len(row)):
                    dict_lists[row[i]] = []  # Initialize empty lists for each column
                count += 1
            elif count == 1 and skipfirst2rows:
                list_types = [row[i] for i in range(len(row))]  # Map index to column name
                count+=1
            else:
                for i in range(len(row)):
                    if list_types[i]=='int':
                        dict_lists[list_categories[i]].append(int(row[i]))  # Append each cell to its corresponding list
                    elif list_types[i]=='float':
                        dict_lists[list_categories[i]].append(float(row[i]))
                    else: # par defaut ca sera str
                        dict_lists[list_categories[i]].append(row[i])
                count += 1
    
    return list_categories,dict_lists,list_types


# %% maintenant on veut implementer ca dans une sorte de routine

def routine_vitesse_phase(f_exc=float, freq_acq=float, general_folder=str, dosavgol=False,W=int,Dt=int,plot=True,index_profile_line=int, xlim_fit = tuple, dcm = float, dpx = float):
    ## load data
    #general_folder = 'F:/Banquise/Vasco/Frigo_pmmh/' +date+ '/'
    #general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
    #general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

    folder_video_demod = general_folder + str(f_exc) +'Hz_'+ str(freq_acq) +'Hz/matData/video_demod_W'+ str(W) +'_Dt' +str(Dt)
    matfile_path = folder_video_demod + '/figdata_complex.mat'
    print('data loading...')
    Data_demod = load_complex_field(matfile_path)
    print('data loaded!')
    complex_field = np.copy(Data_demod.data)
    if plot:
        plt.figure()
        plt.imshow(np.real(complex_field)/np.max(np.abs(complex_field)))
        plt.show()
    ## process
    nb_periods = 4
    nb_pts_per_period = 256
    #index_profile_line = 34 # pour piv du 25/11 c'est 34 (W=32), pour W=64 on met donc 17
    m = nb_periods*nb_pts_per_period
    matrix = np.zeros((m,complex_field.shape[1]))
    for i in range(m):
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))#/np.max(np.abs(complex_field))
        img = np.real(complex_field*np.exp(-2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        if dosavgol==True:
            img = savgol_filter(img,21,3)
        #img = np.real(complex_field*np.exp(2*np.pi*1j*i/nb_pts_per_period))/np.abs(complex_field)
        p = img[index_profile_line,:]
        #plt.imshow(img)
        #plt.pause(0.01)
        matrix[i,:] = p


    # calcul du spatio-temporel du champ démodulé en temps


    X_MAXS = np.arange(matrix.shape[1])
    Y_MAXS = []
    
    if direction==-1:
        matrix = np.flip(matrix,axis=0)
    elif direction==1:
        pass
    print('shape matrix : ',matrix.shape)
    for i in range(matrix.shape[1]):
        indices = find_peaks(matrix[:,i])[0]
        i0 = indices[0]
        i1 = indices[1]
        if len(Y_MAXS)==0:
            Y_MAXS.append(i0)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=nb_pts_per_period/3)):
        #    Y_MAXS.append(i0)
        #elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>nb_pts_per_period/3)):
        #    Y_MAXS.append(i1)
        elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)<=abs(Y_MAXS[-1]-i1))):
            Y_MAXS.append(i0)
        elif ((len(Y_MAXS)!=0)&(abs(Y_MAXS[-1]-i0)>abs(Y_MAXS[-1]-i1))):
            Y_MAXS.append(i1)
        else:
            Y_MAXS.append(np.nan)

            
    ## plot the line to fit on top of the matrix imshow (test)
    if plot:
        plt.figure()
        plt.title('$f_{exc}=$'+str(f_exc)+' Hz')
        plt.imshow(matrix.T,aspect='auto')
        plt.plot(Y_MAXS,np.arange(len(matrix[0,:])),'r')
        plt.colorbar()
        plt.ylabel('position along profile [PIV box units]',fontsize=15)
        plt.xlabel('time [256 px = T = 1/$f_{exc}$]',fontsize=15)
        #plt.ylim(np.nanmin(Y_MAXS),np.nanmax(Y_MAXS))
        plt.ylim(0,55)
        
        plt.show()

    ##  fit the line on the spatio-temporal diagram

    x_px = np.arange(xlim_fit[0],xlim_fit[1],1)
    x_meters = x_px * W/2 * (dcm*1e-2)/dpx

    #x_meters = np.unwrap(x_meters)

    t_px = Y_MAXS[xlim_fit[0]:xlim_fit[1]]
    print(t_px)
    t_sec = np.array(t_px)*(1/f_exc)/nb_pts_per_period


    def linear(x,a,b):
        return a*x+b
    opt,pcov = curve_fit(linear,t_sec,x_meters)
    plt.plot(t_sec,x_meters,'o')
    t_fit = np.linspace(np.min(t_sec),np.max(t_sec),100)
    plt.plot(t_fit,linear(t_fit,opt[0],opt[1]))
    plt.xlabel('Time (seconds)',fontsize=15)
    plt.ylabel('Position of the wave front (meters)',fontsize=15)
    plt.title('Slope given by the fit : $v_{\phi}$ = '+str(np.round(opt[0],2))+' +- '+str(np.round(np.sqrt(pcov[0][0]),2))+' m/s',fontsize=15)
    plt.show()
    print('phase velocity = '+str(opt[0])+' +- '+str(np.sqrt(pcov[0][0]))+' m/s')
    # il faut bien prendre sqrt(pcov) pour prendre l'ecart type et non la variance !!!
    return opt[0],np.sqrt(pcov[0][0])

#%% changer generalfolder si besoin

#date = '20241203'
date = '0506'

acq_num = 2
camera_SN = '40300722'
#camera_SN = '40437120'


#W = 64
#Dt = 20
W = 64
Dt = 50
direction = -1

computer = 'dell_vasco'

if computer=='dell_vasco':
    general_folder = f'D:/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
elif computer=='Leyre':
    general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

# pour gre 25 sur Babasse :
general_folder = f'R:/Gre25/{date}/cameras/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'



#general_folder = f'F:/manip_grenoble2024/manips_relation_dispersion/{date}/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#general_folder = f'G:/Grenoble/{date}/manip_relation_dispersion/Acquisition_1/camera_{camera_SN}/'

# pour la manip plexi :
#camera_SN = '40307970'
#camera_SN = '40300722'
#date = '20250130'
#acq_num = 5
#general_folder = f'/run/user/1003/gvfs/afp-volume:host=thiou.local,user=vasco,volume=labshared1/Banquise/Vasco/plexi/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
#W = 64
#Dt = 1
#direction = 1

#%% lister toutes les frequences
def list_all_freq(general_folder=general_folder):
    ll = os.listdir(general_folder)
    ld = []
    list_f_exc = []
    list_freq_acq = []
    for l in ll:
        if ('Hz_' in l)&(l[-2:]=='Hz'):
            ld.append(l)
            print(l)
            list_f_exc.append(int(l[:l.find('Hz')]))
            if '.' in l[l.find('_')+1:-2]:
                list_freq_acq.append(float(l[l.find('_')+1:-2]))
            elif ('.' in l[l.find('_')+1:-2])==False:
                list_freq_acq.append(int(l[l.find('_')+1:-2]))

    tab_f_exc = np.array(list_f_exc)
    tab_freq_acq = np.array(list_freq_acq,dtype='object')
    indsort = np.argsort(tab_f_exc)
    tab_f_exc = tab_f_exc[indsort]
    tab_freq_acq = tab_freq_acq[indsort]

    return tab_f_exc,tab_freq_acq

tab_f_exc,tab_freq_acq = list_all_freq(general_folder=general_folder)

#%% lancer routine
tab_v_phase = np.zeros(len(tab_f_exc))
tab_v_phase_err = np.zeros(len(tab_f_exc))
for i in range(len(tab_f_exc)):
    tab_v_phase[i],tab_v_phase_err[i] = routine_vitesse_phase(f_exc=tab_f_exc[i],freq_acq=tab_freq_acq[i],general_folder=general_folder,W=W,Dt=Dt,index_profile_line=3,xlim_fit=(0,45),dcm=16,dpx=1452)#,camera_SN='40437120')#'40300722')#)

#%%  relation dispersion
def v_phase_flexural(f,D):
    rho=1e3
    omega = 2*np.pi*f
    return (D/rho)**(1/5) * omega**(3/5)


E = 3e9
h = 3e-3
nu = 0.4
D = (E*h**3)/(12*(1-nu**2))

mask = (tab_v_phase_err < 1/30*tab_v_phase) & (tab_v_phase>0)# & (tab_f_exc<100)

popt,pcov = curve_fit(v_phase_flexural, tab_f_exc[mask],tab_v_phase[mask])

plt.errorbar(tab_f_exc[mask],tab_v_phase[mask],tab_v_phase_err[mask],marker='.',linestyle='')
plt.plot(tab_f_exc,v_phase_flexural(tab_f_exc,D))
plt.plot(tab_f_exc,v_phase_flexural(tab_f_exc,popt[0]),label='fit : D='+str(np.round(popt[0],1))+' +- '+str(np.round(np.sqrt(pcov[0][0])))+' J')
plt.fill_between(tab_f_exc,v_phase_flexural(tab_f_exc,popt[0]-np.sqrt(pcov[0][0])),v_phase_flexural(tab_f_exc,popt[0]+np.sqrt(pcov[0][0])),color='red',alpha=0.2)

#plt.ylim(0,60)
#plt.xlim(0,200)
plt.xlabel('frequency (Hz)',fontsize=15)
plt.ylabel('$v_{\phi} (m/s)$',fontsize=15)
plt.legend()
plt.show()



# %% creation d'un dictionnaire pour ranger les résultats obtenus

# inputs :

list_categories,dict_lists,list_types = load_csv_input_RDD()# le file path est par defaut celui dans hublot

def assign_input_param2global_variables(dict_lists=dict):
    """
    cree les variables dont les noms et les valeurs sont stockées dans dict_lists en tant que variables globales
    """
    for key in dict_lists:
        globals()[key] = dict_lists[key]

assign_input_param2global_variables(dict_lists=dict_lists)
list_xlim_fit = []
for i in range(len(x_inf_fit)):
    list_xlim_fit.append((x_inf_fit[i],x_sup_fit[i]))
# function to run and store the routine for a chosen acquisition and desired values of params
def routine_et_structure(idx_routine):
    date = list_dates[idx_routine]
    acq_num = list_acq_num[idx_routine]
    W = list_W[idx_routine]
    Dt = list_Dt[idx_routine]
    camera_SN = list_camera_SN[idx_routine]
    index_profile_line = list_indices_profile_line[idx_routine]
    xlim_fit = list_xlim_fit[idx_routine]
    dcm = list_dcm[idx_routine]
    dpx = list_dpx[idx_routine]
    

    if (date in dico) == False:
        dico[date] = {}
    if ('acq'+str(acq_num) in dico[date])==False:
        dico[date]['acq'+str(acq_num)] = {}
    dico_acq = dico[date]['acq'+str(acq_num)]
    if ('camera_'+camera_SN in dico_acq)==False:
        dico_acq['camera_'+camera_SN] = {}
    dico_camera = dico_acq['camera_'+camera_SN]
    if (f'W{str(W)}_Dt{str(Dt)}' in dico_camera)==False:
        dico_camera[f'W{str(W)}_Dt{str(Dt)}'] = {}
    dico_piv = dico_camera[f'W{str(W)}_Dt{str(Dt)}']
    if (f'index_profile_line{str(index_profile_line)}' in dico_piv)==False:
        dico_piv[f'index_profile_line{str(index_profile_line)}'] = {}
    dico_piv_profile = dico_piv[f'index_profile_line{str(index_profile_line)}']
    if (f'_xlim_fit{str(xlim_fit)}' in dico_piv_profile)==False:
        dico_piv_profile[f'_xlim_fit{str(xlim_fit)}'] = {}
    
    dico_pivprof_xlimfit = dico_piv_profile[f'_xlim_fit{str(xlim_fit)}']


    # run the routine if data is not in the dictionary:
    if ('tab_v_phase' in dico_pivprof_xlimfit)==False:
        dico_pivprof_xlimfit['dcm'] = list_dcm[idx_routine]
        dico_pivprof_xlimfit['dpx'] = list_dpx[idx_routine]
        if computer=='DellVasco':
            general_folder = f'K:/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'
        elif computer=='Leyre':
            general_folder = f'/run/user/1003/gvfs/smb-share:server=adour.local,share=hublot24/Gre24/Data/{date}/manip_relation_dispersion/Acquisition_{str(acq_num)}/camera_{camera_SN}/'

        tab_f_exc,tab_freq_acq = list_all_freq(general_folder=general_folder)
        dico_pivprof_xlimfit['tab_f_exc'] = tab_f_exc
        dico_pivprof_xlimfit['tab_freq_acq'] = tab_freq_acq
        tab_v_phase = np.zeros(len(tab_f_exc))
        tab_v_phase_err = np.zeros(len(tab_f_exc))
        for i in range(len(tab_f_exc)):
            tab_v_phase[i],tab_v_phase_err[i] = routine_vitesse_phase(f_exc=tab_f_exc[i],freq_acq=tab_freq_acq[i],general_folder=general_folder,W=W,Dt=Dt,index_profile_line=index_profile_line,xlim_fit=xlim_fit,dcm=dcm,dpx=dpx,plot=False)#,camera_SN='40437120')#'40300722')#)

        dico_pivprof_xlimfit['tab_v_phase'] = tab_v_phase
        dico_pivprof_xlimfit['tab_v_phase_err'] = tab_v_phase_err

#############################################################
# run for several acqusitions
## first, if required, create an empty dict
if ('dico' in globals())==False:
    dico = {}
## then fill it with the results
for idx_routine in range(len(list_dates)):
    routine_et_structure(idx_routine)

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

