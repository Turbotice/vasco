import h5py

#general_folder = 'R:/Gre25/Data/0516/cameras/manip_relation_dispersion/Acquisition_3/camera_40300722/70Hz_157.001Hz/'
#path2data = general_folder + 'matData/'
#matfile = f'{path2data}PIV_processed_i0{0}_N{0}_Dt{50}_b1_W{64}_full_total_processed.mat'


def load_piv_matdata(matfile):
    with h5py.File(matfile, 'r') as fmat:
        mat_dict = {}
        
        print('Top-level keys : ', list(fmat.keys()))

        mat_dict = mat_to_dict(fmat['m'],fmat['m'])
    return mat_dict