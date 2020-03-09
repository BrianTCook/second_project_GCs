import gzip
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

import glob

def maps(code_name, orbiter_name, Norbiters):
    
    print('!!!!!!!!!!!!!!!!!!!!')
    print('now mapping: %s, %s'%(code_name, orbiter_name))
    print('!!!!!!!!!!!!!!!!!!!!')
    
    vel_ticks = [ -200, 0, 200 ]
    pos_ticks = [ -1.5, 0, 1.5 ]
    
    phase_files = glob.glob('all_data_%s_%s_Norbiters_*.npy.gz'%(code_name, orbiter_name))
    
    logN_max = 9
    Norbiters_list = [ int(2**i) for i in range(logN_max)]
    
    plt.rc('font', family='serif')
    
    for filename_phase in phase_files:
    
        for Norbiters in Norbiters_list:
            
            test_str = '_' + str(Norbiters) + '.'
            
            print('filename_phase: %s'%(filename_phase))
            print('test_str: %s'%(test_str))
            
            if test_str in filename_phase:
        
                sim_times = np.loadtxt('times_in_Myr_%s_%s_Norbiters_%i.txt'%(code_name, orbiter_name, Norbiters))    
                colors = np.loadtxt(code_name + '_' + orbiter_name + '_colors_Norbiters_' + str(Norbiters) + '.txt')
                
                f_all = gzip.GzipFile('all_data_%s_%s_Norbiters_%s.npy.gz'%(code_name, orbiter_name, str(Norbiters)), 'r')
                phase_space_data = np.load(f_all)
            
                ntimes, nparticles, ndim = phase_space_data.shape
                colors = np.loadtxt(code_name + '_' + orbiter_name + '_colors_Norbiters_' + str(Norbiters) + '.txt')
            
                for i, t in enumerate(sim_times):
                
                    if i == 0 or i == len(sim_times)-1:
                    
                        w_all = pd.DataFrame(phase_space_data[i,:,:], 
                                                 columns=['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz'])
                            
                        x,y,z = w_all['x'].tolist(), w_all['y'].tolist(), w_all['z'].tolist()
                        vx,vy,vz = w_all['vx'].tolist(), w_all['vy'].tolist(), w_all['vz'].tolist()
                    
                        fig, axs = plt.subplots(5, 5)
                        
                        #first column            
                        
                        axs[0, 0].scatter(x, y, s=1, c=colors)
                        axs[0, 0].set(xlabel='x', ylabel='y')
                        axs[0, 0].set_xlim(-2, 2)
                        axs[0, 0].set_ylim(-2, 2)
                        axs[0, 0].set_xticks(pos_ticks)
                        axs[0, 0].set_yticks(pos_ticks)
                        axs[0, 0].tick_params(labelsize='small')
                        
                        axs[1, 0].scatter(x, z, s=1, c=colors)
                        axs[1, 0].set(xlabel='x', ylabel='z')
                        axs[1, 0].set_xlim(-2, 2)
                        axs[1, 0].set_ylim(-2, 2)
                        axs[1, 0].set_xticks(pos_ticks)
                        axs[1, 0].set_yticks(pos_ticks)
                        axs[1, 0].tick_params(labelsize='small')
                        
                        axs[2, 0].scatter(x, vx,s=1, c=colors)
                        axs[2, 0].set(xlabel='x', ylabel='vx')
                        axs[2, 0].set_xlim(-2, 2)
                        axs[2, 0].set_ylim(-400, 400)
                        axs[2, 0].set_xticks(pos_ticks)
                        axs[2, 0].set_yticks(vel_ticks)
                        axs[2, 0].tick_params(labelsize='small')
                        
                        axs[3, 0].scatter(x, vy, s=1, c=colors)
                        axs[3, 0].set(xlabel='z', ylabel='vy')
                        axs[3, 0].set_xlim(-2, 2)
                        axs[3, 0].set_ylim(-400, 400)
                        axs[3, 0].set_xticks(pos_ticks)
                        axs[3, 0].set_yticks(vel_ticks)
                        axs[3, 0].tick_params(labelsize='small')
                        
                        axs[4, 0].scatter(x, vz, s=1, c=colors)
                        axs[4, 0].set(xlabel='x', ylabel='vz')
                        axs[4, 0].set_xlim(-2, 2)
                        axs[4, 0].set_ylim(-400, 400)
                        axs[4, 0].set_xticks(pos_ticks)
                        axs[4, 0].set_yticks(vel_ticks)
                        axs[4, 0].tick_params(labelsize='small')
                        
                        #second column
                        
                        axs[1, 1].scatter(y, z, s=1, c=colors)
                        axs[1, 1].set(xlabel='y', ylabel='z')
                        axs[1, 1].set_xlim(-2, 2)
                        axs[1, 1].set_ylim(-2, 2)
                        axs[1, 1].set_xticks(pos_ticks)
                        axs[1, 1].set_yticks(pos_ticks)
                        axs[1, 1].tick_params(labelsize='small')
                        
                        axs[2, 1].scatter(y, vx, s=1, c=colors)
                        axs[2, 1].set(xlabel='y', ylabel='vx')
                        axs[2, 1].set_xlim(-2, 2)
                        axs[2, 1].set_ylim(-400, 400)
                        axs[2, 1].set_xticks(pos_ticks)
                        axs[2, 1].set_yticks(vel_ticks)
                        axs[2, 1].tick_params(labelsize='small')
                        
                        axs[3, 1].scatter(y, vy, s=1, c=colors)
                        axs[3, 1].set(xlabel='y', ylabel='vy')
                        axs[3, 1].set_xlim(-2, 2)
                        axs[3, 1].set_ylim(-400, 400)
                        axs[3, 1].set_xticks(pos_ticks)
                        axs[3, 1].set_yticks(vel_ticks)
                        axs[3, 1].tick_params(labelsize='small')
                        
                        axs[4, 1].scatter(y, vz, s=1, c=colors)
                        axs[4, 1].set(xlabel='y', ylabel='vz')
                        axs[4, 1].set_xlim(-2, 2)
                        axs[4, 1].set_ylim(-400, 400)
                        axs[4, 1].set_xticks(pos_ticks)
                        axs[4, 1].set_yticks(vel_ticks)
                        axs[4, 1].tick_params(labelsize='small')
                        
                        #third column
                        
                        axs[2, 2].scatter(z, vx, s=1, c=colors)
                        axs[2, 2].set(xlabel='z', ylabel='vx')
                        axs[2, 2].set_xlim(-2, 2)
                        axs[2, 2].set_ylim(-400, 400)
                        axs[2, 2].set_xticks(pos_ticks)
                        axs[2, 2].set_yticks(vel_ticks)
                        axs[2, 2].tick_params(labelsize='small')
                        
                        axs[3, 2].scatter(z, vy, s=1, c=colors)
                        axs[3, 2].set(xlabel='z', ylabel='vy')
                        axs[3, 2].set_xlim(-2, 2)
                        axs[3, 2].set_ylim(-400, 400)
                        axs[3, 2].set_xticks(pos_ticks)
                        axs[3, 2].set_yticks(vel_ticks)
                        axs[3, 2].tick_params(labelsize='small')
                        
                        axs[4, 2].scatter(z, vz, s=1, c=colors)
                        axs[4, 2].set(xlabel='z', ylabel='vz')
                        axs[4, 2].set_xlim(-2, 2)
                        axs[4, 2].set_ylim(-400, 400)
                        axs[4, 2].set_xticks(pos_ticks)
                        axs[4, 2].set_yticks(vel_ticks)
                        axs[4, 2].tick_params(labelsize='small')
                        
                        #fourth column
                        
                        axs[3, 3].scatter(vx, vy, s=1, c=colors)
                        axs[3, 3].set(xlabel='vx', ylabel='vy')
                        axs[3, 3].set_xlim(-400, 400)
                        axs[3, 3].set_ylim(-400, 400)
                        axs[3, 3].set_xticks(vel_ticks)
                        axs[3, 3].set_yticks(vel_ticks)
                        axs[3, 3].tick_params(labelsize='small')
                                        
                        axs[4, 3].scatter(vx, vz, s=1, c=colors)
                        axs[4, 3].set(xlabel='vx', ylabel='vz')
                        axs[4, 3].set_xlim(-400, 400)
                        axs[4, 3].set_ylim(-400, 400)
                        axs[4, 3].set_xticks(vel_ticks)
                        axs[4, 3].set_yticks(vel_ticks)
                        axs[4, 3].tick_params(labelsize='small')
                        
                        #fifth column
                        
                        axs[4, 4].scatter(vy, vz, s=1, c=colors)
                        axs[4, 4].set(xlabel='vy', ylabel='vz')
                        axs[4, 4].set_xlim(-400, 400)
                        axs[4, 4].set_ylim(-400, 400)
                        axs[4, 4].set_xticks(vel_ticks)
                        axs[4, 4].set_yticks(vel_ticks)
                        axs[4, 4].tick_params(labelsize='small')

                        # Hide x labels and tick labels for top plots and y ticks for right plots.
                        for ax in axs.flat:
                            ax.label_outer()
                            
                        fig.suptitle('Time = %.02f Myr'%(t), fontsize=14)
                        plt.subplots_adjust(wspace=0, hspace=0)
                        plt.savefig('phase_space_map_frame_%s_%s_%s_Norbiters_%s.pdf'%(str(i).rjust(5, '0'), code_name, orbiter_name, str(Norbiters)))
                        plt.close()
                    
                        plt.figure()
                        plt.scatter(x, y, s=1, c=colors)
                        plt.xlim(-2.0, 2.0)
                        plt.ylim(-2.0, 2.0)
                        plt.xlabel('x (kpc)', fontsize=12)
                        plt.ylabel('y (kpc)', fontsize=12)
                        plt.title('time = %.02f Myr'%(t), fontsize=16)
                        plt.savefig('snapshot_%s_%s_%s_Norbiters_%s.pdf'%(code_name, orbiter_name, str(i).rjust(5, '0'), str(Norbiters)))
                        plt.close()
            
    return 0
