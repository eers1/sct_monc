#!/usr/bin/env python3.8
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import re
import sys
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl
import os

def split_at_bl(ds, da, cloud_top_heights, num_bool):
    if num_bool==True:
        da_bl = [da[i].where(ds.z<height).mean()*1e-6 for i, height in enumerate(cloud_top_heights)]
        da_ft = [da[i].where(ds.z>height).mean()*1e-6 for i, height in enumerate(cloud_top_heights)]
    else:
        da_bl = [da[i].where(ds.z<height).mean() for i, height in enumerate(cloud_top_heights)]
        da_ft = [da[i].where(ds.z>height).mean() for i, height in enumerate(cloud_top_heights)]
    return da_bl, da_ft

def min_max(da):
    return da.min().item(), da.max().item()

def plot_aero(frame, ax, aero_max, title, xvis, yvis, cbar_label):
    c_obj = frame.plot(y='z',ax=ax, norm=colors.Normalize(vmin=0, vmax=aero_max),ylim=(0,1800))
    ax.set(ylabel='Height (m)', xlabel='y (km)',title=title)
    ax.xaxis.set_visible(xvis)
    ax.yaxis.set_visible(yvis)
    c_obj.colorbar.set_label(cbar_label)
    return None


path = sys.argv[1]  # path to merged file
ds = xr.open_dataset(path)
ds = cfl.ds_fix_dims(ds)

r_path = re.compile(r'carisma/eers/(\S*)/(\S*)/(\D*)(\d*)/(\S*)_merged.nc')
a = r_path.search(path)
project = a.group(1)
design = a.group(2)
key = a.group(3) + a.group(4)
run_name = a.group(5)

cltop_ave = np.mean(ds.cltop.where(ds['cltop']!=0.0), axis=(1,2))
cltop_ind=[int(np.where(ds.z == min(ds.z[np.where(ds.z-i >= 0.)[0]]))[0]) for i in cltop_ave]
free_aero_m_mean = ds.q_aitken_sol_mass + ds.q_accum_sol_mass + ds.q_coarse_sol_mass 
total_aero_m_mean = free_aero_m_mean + ds.q_active_sol_liquid
free_aero_n_mean = ds.q_aitken_sol_number + ds.q_accum_sol_number + ds.q_coarse_sol_number 

free_aero_m_mean_bl, free_aero_m_mean_ft = split_at_bl(ds, free_aero_m_mean, cltop_ave, False)
total_aero_m_mean_bl, total_aero_m_mean_ft = split_at_bl(ds, total_aero_m_mean, cltop_ave, False)
free_aero_num_mean_bl, free_aero_num_mean_ft = split_at_bl(ds, free_aero_n_mean, cltop_ave, True)

ait_m_bl, ait_m_ft = split_at_bl(ds, ds.q_aitken_sol_mass, cltop_ave, False)
accum_m_bl, accum_m_ft = split_at_bl(ds, ds.q_accum_sol_mass, cltop_ave, False)
coarse_m_bl, coarse_m_ft = split_at_bl(ds, ds.q_coarse_sol_mass, cltop_ave, False)
active_m_bl, active_m_ft = split_at_bl(ds, ds.q_active_sol_liquid, cltop_ave, False)

ait_num_bl, ait_num_ft = split_at_bl(ds, ds.q_aitken_sol_number, cltop_ave, True)
accum_num_bl, accum_num_ft = split_at_bl(ds, ds.q_accum_sol_number, cltop_ave, True)
coarse_num_bl, coarse_num_ft = split_at_bl(ds, ds.q_coarse_sol_number, cltop_ave, True)

xsection=128
aer_lim = (0, 2e-9)
figname='animations/{}/{}_mass_qfields_frame_{}.jpg'

aitken_min, aitken_max = min_max(ds.q_aitken_sol_mass)
accum_min, accum_max = min_max(ds.q_accum_sol_mass)
coarse_min, coarse_max = min_max(ds.q_coarse_sol_mass)

for i in range(len(ds.time_coarse)):
    fig = plt.figure(figsize=(22,15))
    gs1 = GridSpec(2,3,left=0.05,bottom=0.4,top=0.95,right=0.95,hspace=0.15,wspace=0.05) 
    a = fig.add_subplot(gs1[0])
    b = fig.add_subplot(gs1[1],sharey=a)
    c = fig.add_subplot(gs1[2],sharey=a)
    d = fig.add_subplot(gs1[3],sharex=a)
    e = fig.add_subplot(gs1[4],sharex=b)
    f = fig.add_subplot(gs1[5],sharex=c)
    gs2 = GridSpec(1,2,left=0.05,bottom=0.05,top=0.35,right=0.8,hspace=0.15,wspace=0.15) 
    g = fig.add_subplot(gs2[0])
    h = fig.add_subplot(gs2[1])
    plt.rc('axes', titlesize=20)

    plot_aero(ds.q_aitken_sol_mass[i,xsection,:,:], a, aitken_max, 'Aitken soluble', False, True, 'mmr (kg kg^-1)')
    plot_aero(ds.q_accum_sol_mass[i,xsection,:,:], b, accum_max, 'Accumulation soluble', False, False, 'mmr (kg kg^-1)')
    plot_aero(ds.q_coarse_sol_mass[i,xsection,:,:], c, coarse_max, 'Coarse soluble', False, False, 'mmr (kg kg^-1)')
    plot_aero(ds.q_active_sol_liquid[i,xsection,:,:], d, 2e-9, 'Active soluble', True, True, 'mmr (kg kg^-1)')
    plot_aero(ds.q_cloud_liquid_mass[i,xsection,:,:], e, 0.0007, 'Cloud liquid', True, False, 'mmr (kg kg^-1)')
    plot_aero(ds.q_rain_mass[i,xsection,:,:], f, 0.000007, 'Rain', True, False, 'mmr (kg kg^-1)')  
        
    g.plot(ds.time_coarse, ait_m_bl, label='Aitken BL mmr')
    g.plot(ds.time_coarse, ait_m_ft, label='Aitken FT mmr')
    g.plot(ds.time_coarse, accum_m_bl, label='Accumulation BL mmr')
    g.plot(ds.time_coarse, accum_m_ft, label='Accumulation FT mmr')
    g.plot(ds.time_coarse, coarse_m_bl, label='Coarse BL mmr')
    g.plot(ds.time_coarse, coarse_m_ft, label='Coarse FT mmr')
    g.plot(ds.time_coarse, free_aero_m_mean_bl, label='Inactive BL mmr')
    g.plot(ds.time_coarse, free_aero_m_mean_ft, label='Inactive FT mmr')
    g.plot(ds.time_coarse, active_m_bl, label='Active BL mmr')
    g.plot(ds.time_coarse, active_m_ft, label='Active FT mmr')
    g.plot(ds.time_coarse, total_aero_m_mean_bl, label='Total BL mmr')
    g.plot(ds.time_coarse, total_aero_m_mean_ft, label='Total FT mmr')

    g.plot((ds.time_coarse[i],ds.time_coarse[i]),(-0.5e-9, 3e-9),c='black')
    
    cfl.add_diurnal(ds, g, (0, 2.1e-9))
    g.set_ylim(0, 2.1e-9)
    g.set_ylabel('Aerosol mean (kg kg^-1)')
    g.legend(loc=(2.2, 0))
    g.set_xlim(0,72)
    g.set_xlabel("Time (hours)")
    
    h.plot(ds.time_coarse, ait_num_bl, label='Aitken BL number')
    h.plot(ds.time_coarse, ait_num_ft, label='Aitken FT number')
    h.plot(ds.time_coarse, accum_num_bl, label='Accumulation BL number')
    h.plot(ds.time_coarse, accum_num_ft, label='Accumulation FT number')
    h.plot(ds.time_coarse, coarse_m_bl, label='Coarse BL number')
    h.plot(ds.time_coarse, coarse_m_ft, label='Coarse FT number')
    h.plot(ds.time_coarse, free_aero_num_mean_bl, label='Inactive BL number')
    h.plot(ds.time_coarse, free_aero_num_mean_ft, label='Inactive FT number')
    
    h.plot((ds.time_coarse[i],ds.time_coarse[i]),(-5,400),c='black')
    
    cfl.add_diurnal(ds, h, (-1, 320))
    h.set_ylim(-1, 320)
    h.set_ylabel('Aerosol mean (cm^-3)')
    h.set_xlim(0,72)
    h.set_xlabel("Time (hours)")

    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(ds.time_coarse)} complete')
    plt.close()
      
figname='animations/{}/{}_num_qfields_frame_{}.jpg'

for i in range(len(ds.time_coarse)):
    fig = plt.figure(figsize=(22,15))
    gs1 = GridSpec(2,3,left=0.05,bottom=0.4,top=0.95,right=0.95,hspace=0.15,wspace=0.05) 
    a = fig.add_subplot(gs1[0])
    b = fig.add_subplot(gs1[1],sharey=a)
    c = fig.add_subplot(gs1[2],sharey=a)
    d = fig.add_subplot(gs1[3],sharex=a)
    e = fig.add_subplot(gs1[4],sharex=b)
    f = fig.add_subplot(gs1[5],sharex=c)
    gs2 = GridSpec(1,2,left=0.05,bottom=0.05,top=0.35,right=0.8,hspace=0.15,wspace=0.15) 
    g = fig.add_subplot(gs2[0])
    h = fig.add_subplot(gs2[1])
    plt.rc('axes', titlesize=20)
    
    plot_aero(ds.q_aitken_sol_number[i,xsection,:,:]*1e-6, a, 200, 'Aitken soluble', False, True, 'Num conc (cm^-3)')
    plot_aero(ds.q_accum_sol_number[i,xsection,:,:]*1e-6, b, 200, 'Aitken soluble', False, False, 'Num conc (cm^-3)')
    plot_aero(ds.q_coarse_sol_number[i,xsection,:,:]*1e-6, c, 200, 'Coarse soluble', False, False, 'Num conc (cm^-3)')
    plot_aero(ds.q_active_sol_liquid[i,xsection,:,:], d, 2e-9, 'Active soluble', True, True, 'mmr (kg kg^-1)')
    plot_aero(ds.q_cloud_liquid_number[i,xsection,:,:]*1e-6, e, 200, 'Cloud liquid', True, False, 'Num conc (cm^-3)')
    plot_aero(ds.q_rain_number[i,xsection,:,:]*1e-6, f, 0.03, 'Rain', True, False, 'Num conc (cm^-3)')
       
    g.plot(ds.time_coarse, ait_m_bl, label='Aitken BL')
    g.plot(ds.time_coarse, ait_m_ft, label='Aitken FT')
    g.plot(ds.time_coarse, accum_m_bl, label='Accumulation BL')
    g.plot(ds.time_coarse, accum_m_ft, label='Accumulation FT')
    g.plot(ds.time_coarse, coarse_m_bl, label='Coarse BL')
    g.plot(ds.time_coarse, coarse_m_ft, label='Coarse FT')
    g.plot(ds.time_coarse, free_aero_m_mean_bl, label='Inactive BL')
    g.plot(ds.time_coarse, free_aero_m_mean_ft, label='Inactive FT')
    g.plot(ds.time_coarse, active_m_bl, label='Active BL')
    g.plot(ds.time_coarse, active_m_ft, label='Active FT')
    g.plot(ds.time_coarse, total_aero_m_mean_bl, label='Total BL')
    g.plot(ds.time_coarse, total_aero_m_mean_ft, label='Total FT')

    g.plot((ds.time_coarse[i],ds.time_coarse[i]),(-0.5e-9, 3e-9),c='black')
    
    cfl.add_diurnal(ds, g, (0, 2.1e-9))
    g.set_ylim(0, 2.1e-9)
    g.set_ylabel('Aerosol mean (kg kg^-1)')
    g.legend(loc=(2.2, 0))
    g.set_xlim(0,72)
    g.set_xlabel("Time (hours)")
    
    h.plot(ds.time_coarse, ait_num_bl, label='Aitken BL number')
    h.plot(ds.time_coarse, ait_num_ft, label='Aitken FT number')
    h.plot(ds.time_coarse, accum_num_bl, label='Accumulation BL number')
    h.plot(ds.time_coarse, accum_num_ft, label='Accumulation FT number')
    h.plot(ds.time_coarse, coarse_m_bl, label='Coarse BL number')
    h.plot(ds.time_coarse, coarse_m_ft, label='Coarse FT number')
    h.plot(ds.time_coarse, free_aero_num_mean_bl, label='Inactive BL number')
    h.plot(ds.time_coarse, free_aero_num_mean_ft, label='Inactive FT number')
    
    h.plot((ds.time_coarse[i],ds.time_coarse[i]),(-5,400),c='black')
    
    cfl.add_diurnal(ds, h, (-1, 320))
    h.set_ylim(-1, 320)
    h.set_ylabel('Aerosol mean (cm^-3)')
    h.set_xlim(0,72)
    h.set_xlabel("Time (hours)")

    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(ds.time_coarse)} complete')
    plt.close()


os.system(f'find /home/users/eers/sct/animations/{run_name} -name "{run_name}_mass_qfields_frame_*.jpg" | sort -V | xargs cat | ffmpeg -f image2pipe -r 1 -vcodec mjpeg -i - -vcodec libx264 /home/users/eers/sct/animations/{run_name}/{run_name}_qfields_mass.mp4')

os.system(f'find /home/users/eers/sct/animations/{run_name} -name "{run_name}_num_qfields_frame_*.jpg" | sort -V | xargs cat | ffmpeg -f image2pipe -r 1 -vcodec mjpeg -i - -vcodec libx264 /home/users/eers/sct/animations/{run_name}/{run_name}_qfields_num.mp4')