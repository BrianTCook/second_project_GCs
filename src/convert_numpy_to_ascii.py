import numpy as np

code_name, orbiter_name, Norbiters = 'tree', 'SingleCluster', 8
data_directory = '/home/brian/Desktop/second_project_gcs/data/'

all_data = np.load(data_directory+'all_data_%s_%s_Norbiter_%s.npy'%(code_name, orbiter_name, str(Norbiters)))
N_timesteps = len(all_data[:,0,0])

for i in range(N_timesteps):
    
    data_to_keep = all_data[i, :, 0:] #gets rid of mass
    np.savetxt('for_enbid_%s_%s_frame_%s_Norbiters_%s.txt'%(code_name, orbiter_name,
                                                            str(i).rjust(5, '0'), str(Norbiters)))