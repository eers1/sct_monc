#!/usr/bin/env python3.8
import re
import sys
import numpy as np
import xarray as xr
from matplotlib.gridspec import GridSpec
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
from copy import copy
import os

sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl

font = {'family':'sans-serif', 'size'   : 15}
rc('font', **font)

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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

def load_csvs(run_name):
    '''
    Load in csvs from lwp mask calculations
    '''
    path='/home/users/eers/sct/lwp_mask_csvs'
    
    loadedArr = np.loadtxt(f'{path}/{run_name}_lwp_mask_2d.csv')
    lwp_mask = loadedArr.reshape(loadedArr.shape[0], loadedArr.shape[1] // 256, 256)
    height = np.loadtxt(f'{path}/{run_name}_inv_height.csv')
    cloud_frac = np.loadtxt(f'{path}/{run_name}_cloud_frac.csv')
    times = np.loadtxt(f'{path}/{run_name}_times.csv')
    cdnc_mean = np.loadtxt(f'{path}/{run_name}_active_aero.csv')
    
    return lwp_mask, height, cloud_frac, times, cdnc_mean

def make_frame(ds, da, frame, i, transect, lwp_mask, cmap, norm, ticks, tick_pos, times, cloud_frac):
    fig = plt.figure(figsize=(21,6))
    gs1 = GridSpec(1,2,left=0.05,bottom=0.1,top=0.95,right=0.68,hspace=0.1,wspace=0.14) 
    ax_cloud_mmr = fig.add_subplot(gs1[0])
    ax_lwp_mask = fig.add_subplot(gs1[1])
    gs2 = GridSpec(2,1,left=0.75,bottom=0.1,top=0.95,right=0.99,hspace=0.1,wspace=0.12) 
    ax_mean_lwp = fig.add_subplot(gs2[0])
    ax_mean_cf = fig.add_subplot(gs2[1])
    
    t=copy(transect)
    
    cloud_mmr=frame[xsection,:,:]*1000
    c_obj1=cloud_mmr.plot(x='y',y='z',vmin=np.min(da.values)*1000, vmax=np.max(da.values)*1000,ax=ax_cloud_mmr)
    #c_obj1=cl_mmr.mean(axis=0).plot(y='z',vmin=cl_min*1000, vmax=cl_max*1000,ax=a)
    c_obj1.colorbar.set_label('mmr ($g~kg^{-1}$)') #,fontsize=font)
    ax_cloud_mmr.set(xlabel='y (km)', ylabel='Height (m)', title='Cloud liquid cross section')
    #rwp = ds.rwp*1000
    #c_obj1=rwp[i].plot(x='y', ax=a, vmin=0, vmax=10)
    #c_obj1.colorbar.set_label('rwp (g m^-2)') #,fontsize=font)
    #a.set(xlabel='y (km)', ylabel='x (km)', title='Rain water path')
    #print("plot imshow")
    
    c_obj2 = ax_lwp_mask.imshow(lwp_mask[i] ,cmap=cmap, norm=norm)
    sp_mask = ds.surface_precip[i].where(ds.surface_precip[i].values!=0)*1000
    #sp_mask.plot(ax=ax_lwp_mask)
    #ax_lwp_mask.imshow(np.flip(np.array(sp_mask.values), axis=0), cmap="Blues", norm=colors.Normalize(vmin=0, vmax=0.0005))
    time_ints = ds.time_coarse[i]-ds.time_coarse[i-1]
    sp_hr = sp_mask/time_ints
    cs = ax_lwp_mask.contourf(sp_hr, hatches=['.', '..','...','....'],
                  colors='none', levels=[0.0001,0.001,0.01,0.1,1])
    
    for collection in cs.collections:
        collection.set_edgecolor('aqua')
        collection.set_linewidth(0.)
        
    leg_artists = []
    for k in ['.', '..','...','....']:
        p = patches.Patch(facecolor="white", edgecolor='aqua', hatch=k)
        leg_artists.append(p)

    ax_lwp_mask.legend(leg_artists, ["1e-4 $mm~hr^{-1}$","1e-3","1e-2","1e-1","1e0"], loc=(1.05,0), frameon=False)
        
    ax_lwp_mask.invert_yaxis() 
    ax_lwp_mask.add_patch(t)
    ax_lwp_mask.set(xlabel='y (km)', ylabel='x (km)', title='LWP top-down view')
    ax_lwp_mask.set_xticklabels(ticks)
    ax_lwp_mask.set_xticks(tick_pos)
    ax_lwp_mask.set_yticklabels(ticks)
    ax_lwp_mask.set_yticks(tick_pos)
    cbar2=fig.colorbar(c_obj2,ax=ax_lwp_mask, shrink=0.77, anchor=(0.0,0.91))
    cbar2.set_label('LWP ($g~m^{-2}$)')#, size='small')
    #cbar2.ax.tick_params(labelsize='small')
    
    x_lims=(0,da.time_coarse.max())
    lwp_mean=ds.LWP_mean*1000
    rwp_mean=ds.RWP_mean*1000
    lwp_mean.plot(ax=ax_mean_lwp, label='lwp')
    rwp_mean.plot(ax=ax_mean_lwp, label='rwp')
    c_y_lims=(-1,lwp_mean.max()+10)
    ax_mean_lwp.plot((ds.time_coarse[i],ds.time_coarse[i]),c_y_lims,c='black')
    ax_mean_lwp.set(xlim=x_lims,ylim=c_y_lims, ylabel='LWP mean ($g~m^{-2}$)')
    #ax_mean_lwp.set_xlim(0,72) 
    #ax_mean_lwp.set_ylim(0,200) 
    ax_mean_lwp.xaxis.set_visible(False)
    ax_mean_lwp.legend()
    cfl.add_diurnal(ds, ax_mean_lwp, c_y_lims)
    
    ax_mean_cf.plot(times,cloud_frac)
    d_y_lims=(0,1.05)
    ax_mean_cf.plot((ds.time_coarse[i],ds.time_coarse[i]),d_y_lims,c='black')
    ax_mean_cf.set(xlim=x_lims,ylim=d_y_lims, ylabel='Cloud Fraction')
    #ax_mean_cf.set_xlim(0,72)
    ax_mean_cf.set_xlabel("Time (hours from 8LT)")
    cfl.add_diurnal(ds, ax_mean_cf, (-1, 1.05))
    
    #plt.rc('xtick', labelsize='x-small')
    #plt.rc('ytick', labelsize='x-small')
    return fig

# Load ds and csvs
path_to_merged_file = sys.argv[1]
ds, project, design, key, run_name = load_ds(path_to_merged_file)
lwp_mask, height, cloud_frac, times, cdnc_mean = load_csvs(run_name)
da=ds.q_cloud_liquid_mass

# Create a Rectangle patch
xsection=128
transect = patches.Rectangle((0, xsection), len(da.x), 0, linewidth=2, edgecolor='lime', facecolor='none')

# Convert ticks to 2.5km spacing
ticks = np.arange(0,da.x.max(),2.5)
tick_pos=[i*256/da.x.max().values for i in ticks]

# Load LWP cmap 
cmap, norm = cfl.load_lwp_cmap()

figname='animations/{}/{}_lwp_frame_{}.jpg'
print('Starting animation frames')
for i, frame in enumerate(da): 
    fig = make_frame(ds, da, frame, i, transect, lwp_mask, cmap, norm, ticks, tick_pos, times, cloud_frac)
    plt.savefig(figname.format(run_name, run_name, i))
    print(f'Frame {i+1}/{len(da)} complete')
    #plt.show()
    plt.close()
    
# Animate frames
os.system(f'find /home/users/eers/sct/animations/{run_name} -name "{run_name}_lwp_frame_*.jpg" | sort -V | xargs cat | ffmpeg -f image2pipe -r 1 -vcodec mjpeg -i - -vcodec libx264 /home/users/eers/sct/animations/{run_name}/{run_name}_lwp.mp4')