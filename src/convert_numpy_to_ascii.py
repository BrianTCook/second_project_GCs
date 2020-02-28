import gzip
import numpy as np

code_names = [ 'tree', 'Nbody' ]
orbiter_names = [ 'SingleStar', 'SingleCluster' ]
Norbiters_list = [ 1 ]

data_directory = '/home/brian/Desktop/second_project_gcs/data/'

for code_name in code_names:
    for orbiter_name in orbiter_names:
        for Norbiters in Norbiters_list:

            f_all = gzip.GzipFile('all_data_%s_%s_Norbiter_%s.npy.gz'%(code_name, orbiter_name, str(Norbiters)))
            all_data = np.load(f_all)
            N_timesteps = len(all_data[:,0,0])
            
            for i in range(N_timesteps):
                
                data_to_keep = all_data[i, :, 0:] #gets rid of mass
                np.savetxt('for_enbid_%s_%s_frame_%s_Norbiters_%s.txt'%(code_name, orbiter_name,
                                                                        str(i).rjust(5, '0'), str(Norbiters)))