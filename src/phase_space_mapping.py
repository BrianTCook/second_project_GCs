import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

def maps(code_name, orbiter_name):
    
    print('!!!!!!!!!!!!!!!!!!!!')
    print('now mapping: %s, %s'%(code_name, orbiter_name))
    print('!!!!!!!!!!!!!!!!!!!!')
    
    filename_time = 'time_data_%s_%s.npy'%(orbiter_name, code_name)
    filename_phase = 'sixD_data_%s_%s.npy'%(orbiter_name, code_name)    
    
    sim_times = np.load(filename_time)
    phase_space_data = np.load(filename_phase)
    ntimes, ndim, nparticles = phase_space_data.shape

    for i, t in enumerate(sim_times):
    
        w_all = pd.DataFrame(np.rot90(phase_space_data[i,:,:]), 
                             columns=['x', 'y', 'z', 'vx', 'vy', 'vz'])
        
        x,y,z = w_all['x'].tolist(), w_all['y'].tolist(), w_all['z'].tolist()
        vx,vy,vz = w_all['vx'].tolist(), w_all['vy'].tolist(), w_all['vz'].tolist()
        
        '''
        fig, axs = plt.subplots(5, 5)
        
        #first column
        
        axs[0, 0].scatter(x, y, s=2, c='k')
        axs[0, 0].set(xlabel='x', ylabel='y')
        axs[0, 0].set_xlim(-5, 5)
        axs[0, 0].set_ylim(-5, 5)
        
        axs[1, 0].scatter(x, z, s=2, c='k')
        axs[1, 0].set(xlabel='x', ylabel='z')
        axs[1, 0].set_xlim(-5, 5)
        axs[1, 0].set_ylim(-5, 5)
        
        axs[2, 0].scatter(x, vx,s=2, c='k')
        axs[2, 0].set(xlabel='x', ylabel='vx')
        axs[2, 0].set_xlim(-5, 5)
        axs[2, 0].set_ylim(-400, 400)
        
        axs[3, 0].scatter(x, vy, s=2, c='k')
        axs[3, 0].set(xlabel='z', ylabel='vy')
        axs[3, 0].set_xlim(-5, 5)
        axs[3, 0].set_ylim(-400, 400)
        
        axs[4, 0].scatter(x, vz, s=2, c='k')
        axs[4, 0].set(xlabel='x', ylabel='vz')
        axs[4, 0].set_xlim(-5, 5)
        axs[4, 0].set_ylim(-400, 400)
        
        #second column
        
        axs[1, 1].scatter(y, z, s=2, c='k')
        axs[1, 1].set(xlabel='y', ylabel='z')
        axs[1, 1].set_xlim(-5, 5)
        axs[1, 1].set_ylim(-5, 5)
        
        axs[2, 1].scatter(y, vx, s=2, c='k')
        axs[2, 1].set(xlabel='y', ylabel='vx')
        axs[2, 1].set_xlim(-5, 5)
        axs[2, 1].set_ylim(-400, 400)
        
        axs[3, 1].scatter(y, vy, s=2, c='k')
        axs[3, 1].set(xlabel='y', ylabel='vy')
        axs[3, 1].set_xlim(-5, 5)
        axs[3, 1].set_ylim(-400, 400)
        
        axs[4, 1].scatter(y, vz, s=2, c='k')
        axs[4, 1].set(xlabel='y', ylabel='vz')
        axs[4, 1].set_xlim(-5, 5)
        axs[4, 1].set_ylim(-400, 400)
        
        #third column
        
        axs[2, 2].scatter(z, vx, s=2, c='k')
        axs[2, 2].set(xlabel='z', ylabel='vx')
        axs[2, 2].set_xlim(-5, 5)
        axs[2, 2].set_ylim(-400, 400)
        
        axs[3, 2].scatter(z, vy, s=2, c='k')
        axs[3, 2].set(xlabel='z', ylabel='vy')
        axs[3, 2].set_xlim(-5, 5)
        axs[3, 2].set_ylim(-400, 400)
        
        axs[4, 2].scatter(z, vz, s=2, c='k')
        axs[4, 2].set(xlabel='z', ylabel='vz')
        axs[4, 2].set_xlim(-5, 5)
        axs[4, 2].set_ylim(-400, 400)
        
        #fourth column
        
        axs[3, 3].scatter(vx, vy, s=2, c='k')
        axs[3, 3].set(xlabel='vx', ylabel='vy')
        axs[3, 3].set_xlim(-400, 400)
        axs[3, 3].set_ylim(-400, 400)
        
        axs[4, 3].scatter(vx, vz, s=2, c='k')
        axs[4, 3].set(xlabel='vx', ylabel='vz')
        axs[4, 3].set_xlim(-400, 400)
        axs[4, 3].set_ylim(-400, 400)
        
        #fifth column
        
        axs[4, 4].scatter(vy, vz, s=2, c='k')
        axs[4, 4].set(xlabel='vy', ylabel='vz')
        axs[4, 4].set_xlim(-400, 400)
        axs[4, 4].set_ylim(-400, 400)
        
        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()
            
        fig.suptitle('Time = %.02f Myr'%(t), fontsize=14)
        plt.savefig('phase_space_map_frame=%s_%s_%s.png'%(str(i).rjust(4, '0'), code_name, orbiter_name))
        '''
        
        plt.figure()
        plt.scatter(x, y, 'k')
        plt.xlabel('x (kpc)', fontsize=12)
        plt.ylabel('y (kpc)', fontsize=12)
        plt.title('time = %.02f Myr'%(t), fontsize=16)
        plt.savefig('snapshot_%s_%s_%i.png'%(orbiter_name, code_name, i))
        plt.close()
    
    return 0
