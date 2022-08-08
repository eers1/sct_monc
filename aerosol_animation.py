#!/usr/bin/env python3.8
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import re
import sys
from matplotlib.gridspec import GridSpec
from matplotlib import rc
import matplotlib.colors as colors
sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl
import os
import datetime

font = {'family':'sans-serif', 'size'   : 15}
rc('font', **font)

def split_at_bl(ds, da, cloud_top_heights, num_bool):
    if num_bool==True:
        da_bl = [da[i].where(ds.z<height).mean()*1e-6 for i, height in enumerate(cloud_top_heights)]
        da_ft = [da[i].where(ds.z>height).mean()*1e-6 for i, height in enumerate(cloud_top_heights)]
    else:
        da_bl = [da[i].where(ds.z<height).mean() for i, height in enumerate(cloud_top_heights)]
        da_ft = [da[i].where(ds.z>height).mean() for i, height in enumerate(cloud_top_heights)]
    return da_bl, da_ft

def plot_aero(frame, ax, aero_max, xvis, yvis, title=None, ylabel=None, xlabel=None, cbar_label=None):
    c_obj = frame.plot(y='z',ax=ax, norm=colors.Normalize(vmin=0, vmax=aero_max), ylim=(0,1800), cbar_kwargs={'pad':0.02})
    
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel('Height (m)')
        
    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel('y (km)')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title('')

    ax.xaxis.set_visible(xvis)
    ax.yaxis.set_visible(yvis)
    if cbar_label:
        c_obj.colorbar.set_label(cbar_label)
    else:
        c_obj.colorbar.set_label('')
    return None

def load_ds(path):
    '''
    Load a ds from specified path and get details of simulation from path
    '''
    ds = xr.open_dataset(path,engine='netcdf4')
    ds = cfl.ds_fix_dims(ds)
    
    r_path = re.compile(r'carisma/eers/(\S*)/(\S*)/(\D*)(\d*)/(\S*)_merged.nc')
    a = r_path.search(path)
    project = a.group(1)
    design = a.group(2)
    key = a.group(3) + a.group(4)
    run_name = a.group(5)
    return ds, project, design, key, run_name

def load_aero_csvs(run_name, aero_keys, aero_colours, aero_style):
    '''
    Load in csvs from splitting the atmosphere at the boundary layer
    '''  
    aero_dict = {}
    for key, colour, style in zip(aero_keys, aero_colours, aero_style):
        aero_dict[key] = [np.loadtxt(f'/home/users/eers/sct/aero_ani_csvs/{run_name}_{key}.csv', delimiter=','), colour, style]
        
    return aero_dict

def fontsizes():
    SMALL_SIZE = 12
    MEDIUM_SIZE = 14
    BIGGER_SIZE = 15

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    
    return None

def create_frame(ds, xsection, aero_m_dict, aero_num_dict, labels, figname):
    fig = plt.figure(figsize=(25,11))
    gs1 = GridSpec(3,5,left=0.03,bottom=0.07,top=0.95,right=0.98,hspace=0.13,wspace=0.04) 
    ax_ait_m = fig.add_subplot(gs1[0])
    ax_acc_m = fig.add_subplot(gs1[1])
    ax_coa_m = fig.add_subplot(gs1[4])
    ax_cl_m = fig.add_subplot(gs1[2])
    ax_r_m = fig.add_subplot(gs1[3])
    
    ax_ait_n = fig.add_subplot(gs1[5])
    ax_acc_n = fig.add_subplot(gs1[6])
    ax_coa_n = fig.add_subplot(gs1[9])
    ax_cl_n = fig.add_subplot(gs1[7])
    ax_r_n = fig.add_subplot(gs1[8])
    
    ax_act_m = fig.add_subplot(gs1[10])
    
    gs2 = GridSpec(1,2,left=0.25,bottom=0.07,top=0.35,right=0.87, wspace=0.13) 
    ax_mmr = fig.add_subplot(gs2[0])
    ax_num = fig.add_subplot(gs2[1])

    plot_aero(ds.q_aitken_sol_mass[i,xsection,:,:], ax_ait_m, 2.5e-11, False, True, title='Aitken', ylabel='')
    plot_aero(ds.q_accum_sol_mass[i,xsection,:,:], ax_acc_m, 1e-8, False, False, title='Accumulation')
    plot_aero(ds.q_coarse_sol_mass[i,xsection,:,:], ax_coa_m, 1e-8, False, False, title='Coarse', cbar_label='mmr (kg kg^-1)')
    plot_aero(ds.q_active_sol_liquid[i,xsection,:,:], ax_act_m, 1e-8, True, True, title='Active', ylabel='')
    plot_aero(ds.q_cloud_liquid_mass[i,xsection,:,:], ax_cl_m, 5e-4, False, False, title='Cloud')   # cloud-averaged cloud mmr is ~1e-4
    plot_aero(ds.q_rain_mass[i,xsection,:,:], ax_r_m, 1e-5, False, False, title='Rain')  # rain-averaged rain mmr is ~1e-6
    
    plot_aero(ds.q_aitken_sol_number[i,xsection,:,:]*1e-6, ax_ait_n, 200, False, True, title='')
    plot_aero(ds.q_accum_sol_number[i,xsection,:,:]*1e-6, ax_acc_n, 400, False, False, title='')
    plot_aero(ds.q_coarse_sol_number[i,xsection,:,:]*1e-6, ax_coa_n, 1, False, False, title='', cbar_label='Num conc (cm^-3)')
    plot_aero(ds.q_cloud_liquid_number[i,xsection,:,:]*1e-6, ax_cl_n, 400, False, False,  title='')
    plot_aero(ds.q_rain_number[i,xsection,:,:]*1e-6, ax_r_n, 0.03, False, False, title='')
        
    lines = []
    for val in aero_m_dict.values():
        line, = ax_mmr.plot(ds.time_coarse, val[0], c=val[1],linestyle=val[2])
        lines.append(line)
        
    ylim_mmr = (0,2.5e-8) 
    ax_mmr.plot((ds.time_coarse[i],ds.time_coarse[i]),ylim_mmr,c='black')
    cfl.add_diurnal(ds, ax_mmr, ylim_mmr)
    ax_mmr.set_ylim(ylim_mmr)
    ax_mmr.set_ylabel('Aerosol mean (kg kg^-1)')
    plain = mlines.Line2D([], [], color='black', linestyle='-')
    dotted = mlines.Line2D([], [], color='black', linestyle=':')
    labels.extend(['Boundary layer', 'Free tropo'])
    ax_mmr.legend([lines[0],lines[2],lines[4],lines[6],lines[8],lines[10], plain, dotted], labels, loc=(2.15, 0))
    ax_mmr.set_xlim(0,72)
    ax_mmr.set_xlabel("Time (hours)")
    
    for val in aero_num_dict.values():
        ax_num.plot(ds.time_coarse, val[0], c=val[1],linestyle=val[2])
    
    ylim_num = (-1, 500)
    ax_num.plot((ds.time_coarse[i],ds.time_coarse[i]), ylim_num,c='black')
    cfl.add_diurnal(ds, ax_num, ylim_num)
    ax_num.set_ylim(ylim_num)
    ax_num.set_ylabel('Aerosol mean (cm^-3)')
    ax_num.set_xlim(0,72)
    ax_num.set_xlabel("Time (hours)")

    fontsizes()
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(ds.time_coarse)} complete')
    plt.close()
    
    
# Load ds and csvs
print(f'Start time: {datetime.datetime.now()}')
path_to_merged_file = sys.argv[1]
ds, project, design, key, run_name = load_ds(path_to_merged_file)

print(f"Project: {project}")
print(f"Simulation: {key}")

aero_m_dict = load_aero_csvs(run_name, ['ait_m_bl', 'ait_m_ft', 'accum_m_bl', 'accum_m_ft', 'coarse_m_bl', 'coarse_m_ft', 'free_aero_m_mean_bl', 'free_aero_m_mean_ft', 
             'active_m_bl', 'active_m_ft', 'total_aero_m_mean_bl', 'total_aero_m_mean_ft'], ['C0', 'C0', 'C1', 'C1', 'C2', 'C2', 'C3', 'C3', 'C4', 'C4', 'C5', 'C5'],['-',':']*6)

aero_num_dict = load_aero_csvs(run_name, ['ait_num_bl', 'ait_num_ft', 'accum_num_bl', 'accum_num_ft', 'coarse_num_bl', 'coarse_num_ft', 'free_aero_num_mean_bl', 'free_aero_num_mean_ft'], 
                             ['C0', 'C0', 'C1', 'C1', 'C2', 'C2', 'C3', 'C3'],['-',':']*4)

labels = ['Aitken', 'Accumulation', 'Coarse', 'Inactive', 'Active', 'Total']

# Create frames 
xsection = 128
figname ='animations/{}/{}_qfields_frame_{}.jpg'
print('Starting animation frames')
for i in range(len(ds.time_coarse)):
    create_frame(ds, xsection, aero_m_dict, aero_num_dict, labels, figname)

# Animate frames
print("Animating")
os.system(f'find /home/users/eers/sct/animations/{run_name} -name "{run_name}_qfields_frame_*.jpg" | sort -V | xargs cat | ffmpeg -f image2pipe -r 1 -vcodec mjpeg -i - -vcodec libx264 /home/users/eers/sct/animations/{run_name}/{run_name}_qfields.mp4')

print(f'End time: {datetime.datetime.now()}')