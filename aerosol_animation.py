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

free_aero_m_mean_bl=[free_aero_m_mean[i].where(ds.z<height).mean() for i, height in enumerate(cltop_ave)]
free_aero_m_mean_ft=[free_aero_m_mean[i].where(ds.z>height).mean() for i, height in enumerate(cltop_ave)]
total_aero_m_mean_bl=[total_aero_m_mean[i].where(ds.z<height).mean() for i, height in enumerate(cltop_ave)]
total_aero_m_mean_ft=[total_aero_m_mean[i].where(ds.z>height).mean() for i, height in enumerate(cltop_ave)]
free_aero_num_mean_bl=[free_aero_n_mean[i].where(ds.z<height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
free_aero_num_mean_ft=[free_aero_n_mean[i].where(ds.z>height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
    
ait_m_bl=[ds.q_aitken_sol_mass[i].where(ds.z<height).mean() for i, height in enumerate(cltop_ave)]
ait_m_ft=[ds.q_aitken_sol_mass[i].where(ds.z>height).mean() for i, height in enumerate(cltop_ave)]
accum_m_bl=[ds.q_accum_sol_mass[i].where(ds.z<height).mean() for i, height in enumerate(cltop_ave)]
accum_m_ft=[ds.q_accum_sol_mass[i].where(ds.z>height).mean() for i, height in enumerate(cltop_ave)]
coarse_m_bl=[ds.q_coarse_sol_mass[i].where(ds.z<height).mean() for i, height in enumerate(cltop_ave)]
coarse_m_ft=[ds.q_coarse_sol_mass[i].where(ds.z>height).mean() for i, height in enumerate(cltop_ave)]
active_m_bl=[ds.q_active_sol_liquid[i].where(ds.z<height).mean() for i, height in enumerate(cltop_ave)]
active_m_ft=[ds.q_active_sol_liquid[i].where(ds.z>height).mean() for i, height in enumerate(cltop_ave)]

ait_num_bl=[ds.q_aitken_sol_number[i].where(ds.z<height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
ait_num_ft=[ds.q_aitken_sol_number[i].where(ds.z>height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
accum_num_bl=[ds.q_accum_sol_number[i].where(ds.z<height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
accum_num_ft=[ds.q_accum_sol_number[i].where(ds.z>height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
coarse_num_bl=[ds.q_coarse_sol_number[i].where(ds.z<height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
coarse_num_ft=[ds.q_coarse_sol_number[i].where(ds.z>height).mean()*1e-6 for i, height in enumerate(cltop_ave)]
xsection=128
aer_lim = (0, 2e-9)
figname='animations/{}/{}_mass_qfields_frame_{}.jpg'

aitken_min=ds.q_aitken_sol_mass.min().item()
aitken_max=ds.q_aitken_sol_mass.max().item()
accum_min=ds.q_accum_sol_mass.min().item()
accum_max=ds.q_accum_sol_mass.max().item()
coarse_min=ds.q_coarse_sol_mass.min().item()
coarse_max=ds.q_coarse_sol_mass.max().item()

for i, frame in enumerate(ds.q_cloud_liquid_mass):
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

    aitken_col = ds.q_aitken_sol_mass[i,xsection,:,:].plot(y='z',ax=a, norm=colors.Normalize(vmin=0, vmax=aitken_max),ylim=(0,1800))
    a.set(ylabel='Height (m)', xlabel='y (km)',title='Aitken soluble')
    a.xaxis.set_visible(False)
    aitken_col.colorbar.set_label('mmr (kg kg^-1)') 
    
    accum_col = ds.q_accum_sol_mass[i,xsection,:,:].plot(y='z',ax=b, norm=colors.Normalize(vmin=0, vmax=accum_max),ylim=(0,1800))
    b.set(ylabel='Height (m)', xlabel='y (km)',title='Accumulation soluble')
    b.xaxis.set_visible(False)
    b.yaxis.set_visible(False)
    accum_col.colorbar.set_label('mmr (kg kg^-1)') 
    
    coarse_col = ds.q_coarse_sol_mass[i,xsection,:,:].plot(y='z',ax=c, norm=colors.Normalize(vmin=0, vmax=coarse_max),ylim=(0,1800))
    c.set(ylabel='Height (m)', xlabel='y (km)',title='Coarse soluble')
    c.xaxis.set_visible(False)
    c.yaxis.set_visible(False)
    coarse_col.colorbar.set_label('mmr (kg kg^-1)') 
    
    active_col = ds.q_active_sol_liquid[i,xsection,:,:].plot(y='z',ax=d, norm=colors.Normalize(vmin=0, vmax=2e-9),ylim=(0,1800))
    d.set(ylabel='Height (m)', xlabel='y (km)',title='Active soluble ')
    active_col.colorbar.set_label('mmr (kg kg^-1)') 
    
    cloud_col = ds.q_cloud_liquid_mass[i,xsection,:,:].plot(y='z',ax=e, norm=colors.Normalize(vmin=0, vmax=0.0007),ylim=(0,1800))
    e.set(ylabel='Height (m)', xlabel='y (km)',title='Cloud liquid')
    e.yaxis.set_visible(False)
    cloud_col.colorbar.set_label('mmr (kg kg^-1)') 
    
    rain_col = ds.q_rain_mass[i,xsection,:,:].plot(y='z',ax=f, norm=colors.Normalize(vmin=0, vmax=0.000007),ylim=(0,1800))
    f.set(ylabel='Height (m)', xlabel='y (km)',title='Rain')
    rain_col.colorbar.set_label('mmr (kg kg^-1)') 
    f.yaxis.set_visible(False)
    
#     g.plot(ds.time_coarse, ts_dict['ait_m_bl'], label='Aitken BL mmr')
#     g.plot(ds.time_coarse, ts_dict['ait_m_ft'], label='Aitken FT mmr')
#     g.plot(ds.time_coarse, ts_dict['accum_m_bl'], label='Accumulation BL mmr')
#     g.plot(ds.time_coarse, ts_dict['accum_m_ft'], label='Accumulation FT mmr')
#     g.plot(ds.time_coarse, ts_dict['coarse_m_bl'], label='Coarse BL mmr')
#     g.plot(ds.time_coarse, ts_dict['coarse_m_ft'], label='Coarse FT mmr')
#     g.plot(ds.time_coarse, ts_dict['active_m_bl'], label='Active BL mmr')
#     g.plot(ds.time_coarse, ts_dict['active_m_ft'], label='Active FT mmr')
#     g.plot(ds.time_coarse, ts_dict['free_aero_m_mean_bl'], label='Inactive BL mmr')
#     g.plot(ds.time_coarse, ts_dict['free_aero_m_mean_ft'], label='Inactive FT mmr')
#     g.plot(ds.time_coarse, ts_dict['total_aero_m_mean_bl'], label='Total BL mmr')
#     g.plot(ds.time_coarse, ts_dict['total_aero_m_mean_ft'], label='Total FT mmr')
    
        
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

    g.plot((ds.time_coarse[i],ds.time_coarse[i]),(-0.5e-9, 3e-9),c='lime')
    
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
    
#     h.plot(ds.time_coarse, ts_dict['ait_num_bl'], label='Aitken BL number')
#     h.plot(ds.time_coarse, ts_dict['ait_num_ft'], label='Aitken FT number')
#     h.plot(ds.time_coarse, ts_dict['accum_num_bl'], label='Accumulation BL number')
#     h.plot(ds.time_coarse, ts_dict['accum_num_ft'], label='Accumulation FT number')
#     h.plot(ds.time_coarse, ts_dict['coarse_num_bl'], label='Coarse BL number')
#     h.plot(ds.time_coarse, ts_dict['coarse_num_ft'], label='Coarse FT number')
#     h.plot(ds.time_coarse, ts_dict['free_aero_num_mean_bl'], label='Inactive BL number')
#     h.plot(ds.time_coarse, ts_dict['free_aero_num_mean_ft'], label='Inactive FT number')
    
    h.plot((ds.time_coarse[i],ds.time_coarse[i]),(-5,400),c='lime')
    
    cfl.add_diurnal(ds, h, (-1, 320))
    h.set_ylim(-1, 320)
    h.set_ylabel('Aerosol mean (cm^-3)')
    h.set_xlim(0,72)
    h.set_xlabel("Time (hours)")
    #h.legend(loc=)

    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(da)} complete')
    #plt.show()
    plt.close()
    
    
figname='animations/{}/{}_num_qfields_frame_{}.jpg'
# aitken_min=ds.q_aitken_sol_number.min().item()*1e-6
# aitken_max=ds.q_aitken_sol_number.max().item()*1e-6
# accum_min=ds.q_accum_sol_number.min().item()*1e-6
# accum_max=ds.q_accum_sol_number.max().item()*1e-6
# coarse_min=ds.q_coarse_sol_number.min().item()*1e-6
# coarse_max=ds.q_coarse_sol_number.max().item()*1e-6
#active_min=ds.q_active_sol_liquid.min().item()
#active_max=ds.q_active_sol_liquid.max().item()
#cloud_min=ds.q_cloud_liquid_mass.min().item()
#cloud_max=ds.q_cloud_liquid_mass.max().item()
#rain_min=ds.q_rain_mass.min().item()
#rain_max=ds.q_rain_mass.max().item()

for i, frame in enumerate(ds.q_cloud_liquid_number):
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

    aitken = ds.q_aitken_sol_number[i,xsection,:,:]*1e-6
    aitken_col=aitken.plot(y='z',ax=a, norm=colors.Normalize(vmin=0, vmax=200),ylim=(0,1800))
    a.set(ylabel='Height (m)', xlabel='y (km)',title='Aitken soluble')
    a.xaxis.set_visible(False)
    aitken_col.colorbar.set_label('Num conc (cm^-3)') 
    
    accum = ds.q_accum_sol_number[i,xsection,:,:]*1e-6
    accum_col = accum.plot(y='z',ax=b, norm=colors.Normalize(vmin=0, vmax=200),ylim=(0,1800))
    b.set(ylabel='Height (m)', xlabel='y (km)',title='Accumulation soluble')
    b.xaxis.set_visible(False)
    b.yaxis.set_visible(False)
    accum_col.colorbar.set_label('Num conc (cm^-3)') 
    
    coarse = ds.q_coarse_sol_number[i,xsection,:,:]*1e-6
    coarse_col = coarse.plot(y='z',ax=c, norm=colors.Normalize(vmin=0, vmax=0.03),ylim=(0,1800))
    c.set(ylabel='Height (m)', xlabel='y (km)',title='Coarse soluble')
    c.xaxis.set_visible(False)
    c.yaxis.set_visible(False)
    coarse_col.colorbar.set_label('Num conc (cm^-3)') 
    
    active = ds.q_active_sol_liquid[i,xsection,:,:]
    active_col = active.plot(y='z',ax=d, norm=colors.Normalize(vmin=0, vmax=2e-9),ylim=(0,1800))
    d.set(ylabel='Height (m)', xlabel='y (km)',title='Active soluble ')
    active_col.colorbar.set_label('mmr (kg kg^-1)') 
    
    cloud = ds.q_cloud_liquid_number[i,xsection,:,:]*1e-6
    cloud_col = cloud.plot(y='z',ax=e, norm=colors.Normalize(vmin=0, vmax=200),ylim=(0,1800))
    e.set(ylabel='Height (m)', xlabel='y (km)',title='Cloud liquid')
    e.yaxis.set_visible(False)
    cloud_col.colorbar.set_label('Num conc (cm^-3)') 
    
    rain = ds.q_rain_number[i,xsection,:,:]*1e-6
    rain_col = rain.plot(y='z',ax=f, norm=colors.Normalize(vmin=0, vmax=0.03),ylim=(0,1800))
    f.set(ylabel='Height (m)', xlabel='y (km)',title='Rain')
    rain_col.colorbar.set_label('Num conc (cm^-3)') 
    f.yaxis.set_visible(False)
    
#     g.plot(ds.time_coarse, ts_dict['ait_m_bl'], label='Aitken BL mmr')
#     g.plot(ds.time_coarse, ts_dict['ait_m_ft'], label='Aitken FT mmr')
#     g.plot(ds.time_coarse, ts_dict['accum_m_bl'], label='Accumulation BL mmr')
#     g.plot(ds.time_coarse, ts_dict['accum_m_ft'], label='Accumulation FT mmr')
#     g.plot(ds.time_coarse, ts_dict['coarse_m_bl'], label='Coarse BL mmr')
#     g.plot(ds.time_coarse, ts_dict['coarse_m_ft'], label='Coarse FT mmr')
#     g.plot(ds.time_coarse, ts_dict['active_m_bl'], label='Active BL mmr')
#     g.plot(ds.time_coarse, ts_dict['active_m_ft'], label='Active FT mmr')
#     g.plot(ds.time_coarse, ts_dict['free_aero_m_mean_bl'], label='Inactive BL mmr')
#     g.plot(ds.time_coarse, ts_dict['free_aero_m_mean_ft'], label='Inactive FT mmr')
#     g.plot(ds.time_coarse, ts_dict['total_aero_m_mean_bl'], label='Total BL mmr')
#     g.plot(ds.time_coarse, ts_dict['total_aero_m_mean_ft'], label='Total FT mmr')
    
        
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

    g.plot((ds.time_coarse[i],ds.time_coarse[i]),(-0.5e-9, 3e-9),c='lime')
    
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
    
#     h.plot(ds.time_coarse, ts_dict['ait_num_bl'], label='Aitken BL number')
#     h.plot(ds.time_coarse, ts_dict['ait_num_ft'], label='Aitken FT number')
#     h.plot(ds.time_coarse, ts_dict['accum_num_bl'], label='Accumulation BL number')
#     h.plot(ds.time_coarse, ts_dict['accum_num_ft'], label='Accumulation FT number')
#     h.plot(ds.time_coarse, ts_dict['coarse_num_bl'], label='Coarse BL number')
#     h.plot(ds.time_coarse, ts_dict['coarse_num_ft'], label='Coarse FT number')
#     h.plot(ds.time_coarse, ts_dict['free_aero_num_mean_bl'], label='Inactive BL number')
#     h.plot(ds.time_coarse, ts_dict['free_aero_num_mean_ft'], label='Inactive FT number')
    
    h.plot((ds.time_coarse[i],ds.time_coarse[i]),(-5,400),c='lime')
    
    cfl.add_diurnal(ds, h, (-1, 320))
    h.set_ylim(-1, 320)
    h.set_ylabel('Aerosol mean (cm^-3)')
    h.set_xlim(0,72)
    h.set_xlabel("Time (hours)")
    #h.legend(loc=)

    plt.rc('xtick', labelsize='small')
    plt.rc('ytick', labelsize='small')
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(da)} complete')
    #plt.show()
    plt.close()
