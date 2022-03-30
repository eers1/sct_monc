#!/usr/bin/env python3.8
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches
import matplotlib.colors as colors
import re
import sys
from copy import copy
sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl

path = sys.argv[1]  # path to merged file
ds = xr.open_dataset(path,engine='netcdf4')
ds = cfl.ds_fix_dims(ds)

r_path = re.compile(r'carisma/eers/(\S*)/(\S*)/(\D*)(\d*)/(\S*)_merged.nc')
a = r_path.search(path)
project = a.group(1)
design = a.group(2)
key = a.group(3) + a.group(4)
run_name = a.group(5)

path='/home/users/eers/sct/lwp_mask_csvs'

loadedArr = np.loadtxt(f'{path}/{run_name}_lwp_mask_2d.csv')
lwp_mask = loadedArr.reshape(loadedArr.shape[0], loadedArr.shape[1] // 256, 256)
height = np.loadtxt(f'{path}/{run_name}_inv_height.csv')
cloud_frac = np.loadtxt(f'{path}/{run_name}_cloud_frac.csv')
times = np.loadtxt(f'{path}/{run_name}_times.csv')
cdnc_mean = np.loadtxt(f'{path}/{run_name}_active_aero.csv')
#arr_dict[name] = [ds, loadedOriginal,inv_height,cloud_frac,cf_times,active_aero]

print('Starting animation frames')
#ds=DS['base']
#key='sct_em0'
#arrays=arr_dict['sct_em0']
figname='animations/{}/{}_frame_{}.jpg'
da=ds.q_cloud_liquid_mass # [12:40]
#print("finding min max")
cl_max=np.max(da.values)
cl_min=np.min(da.values)
rwp_max=10
rwp_min=0
xsection=128

#Create a Rectangle patch
transect = patches.Rectangle((0, xsection), len(da.x), 0, linewidth=2, edgecolor='lime', facecolor='none')
ticks = np.arange(0,da.x.max(),2.5)
tick_pos=[i*256/da.x.max().values for i in ticks]
font=15

cmap, norm = cfl.load_lwp_cmap()

for i, frame in enumerate(da): 
    fig = plt.figure(figsize=(21,6))
    gs1 = GridSpec(1,2,left=0.05,bottom=0.1,top=0.9,right=0.7,hspace=0.1,wspace=0.1) 
    a = fig.add_subplot(gs1[0])
    b = fig.add_subplot(gs1[1])
    gs2 = GridSpec(2,1,left=0.75,bottom=0.1,top=0.9,right=0.95,hspace=0.1,wspace=0.12) 
    c = fig.add_subplot(gs2[0])
    d = fig.add_subplot(gs2[1])
    
    t=copy(transect)
    cl_mmr=frame[xsection,:,:]*1000
    #cl_mmr=frame[:,:,:]*1000
    #print("plotting cross section")
    c_obj1=cl_mmr.plot(x='y',y='z',vmin=cl_min*1000, vmax=cl_max*1000,ax=a)
    #c_obj1=cl_mmr.mean(axis=0).plot(y='z',vmin=cl_min*1000, vmax=cl_max*1000,ax=a)
    c_obj1.colorbar.set_label('mmr ($g~kg^{-1}$)') #,fontsize=font)
    a.set(xlabel='y (km)', ylabel='Height (m)', title='Cloud liquid cross section')
    #rwp = ds.rwp*1000
    #c_obj1=rwp[i].plot(x='y', ax=a, vmin=rwp_min, vmax=rwp_max)
    #c_obj1.colorbar.set_label('rwp (g m^-2)') #,fontsize=font)
    #a.set(xlabel='y (km)', ylabel='x (km)', title='Rain water path')
    #print("plot imshow")
    c_obj2=b.imshow(lwp_mask[i] ,cmap=cmap, norm=norm)
    b.invert_yaxis() 
    b.add_patch(t)
    b.set(xlabel='y (km)', ylabel='x (km)', title='LWP top-down view')
    b.set_xticklabels(ticks)
    b.set_xticks(tick_pos)
    b.set_yticklabels(ticks)
    b.set_yticks(tick_pos)
    cbar2=fig.colorbar(c_obj2,ax=b)
    cbar2.set_label('LWP ($g~m^{-2}$)')
    #print("plotting WPs")
    x_lims=(0,da.time_coarse.max())
    #c.plot(times, total_lwp)
    #lwp_mean=ds.LWP_mean[100:180]*1000
    #rwp_mean=ds.RWP_mean[100:180]*1000
    lwp_mean=ds.LWP_mean*1000
    rwp_mean=ds.RWP_mean*1000
    lwp_mean.plot(ax=c, label='lwp')
    rwp_mean.plot(ax=c, label='rwp')
    c_y_lims=(-1,lwp_mean.max()+10)
    c.plot((da.time_coarse[i],da.time_coarse[i]),c_y_lims,c='lime')
    c.set(xlim=x_lims,ylim=c_y_lims, ylabel='LWP mean ($g~m^{-2}$)')
    c.set_xlim(0,72) # c.set_xlim(0,72)
    c.xaxis.set_visible(False)
    c.legend()
    cfl.add_diurnal(ds, c, c_y_lims)
    #print("plotting cloud frac")
    #d.plot(times[12:40], cloud_frac[12:40])
    d.plot(times,cloud_frac)
    d_y_lims=(0,1.05)
    d.plot((da.time_coarse[i],da.time_coarse[i]),d_y_lims,c='lime')
    d.set(xlim=x_lims,ylim=d_y_lims, ylabel='Cloud Fraction')
    #d.set_xlim(10,25) 
    d.set_xlim(0,72)
    d.set_xlabel("Time (hours from 8LT)")
    cfl.add_diurnal(ds, d, (-1, 1.05))
    
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(da)} complete')
    #plt.show()
    plt.close()
