import numpy as np


def correction_factor(v,ypix,ypix_center,dcm):
    dfactor_sur_dcm = 8e-4
    dcm_sur_dpx = 0.07 # en gros
    dfactor_sur_dpx = dfactor_sur_dcm/dcm_sur_dpx

    dy = ypix-ypix_center
    

    return v * (1/np.cos(alpha_rad))
