#!/usr/bin/env python3.8
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors

def ds_fix_dims(ds):
    ds = ds.rename({ds[ds.w.dims[0]].name:'time_coarse', ds[ds.w_max.dims[0]].name:'time_fine', ds[ds.rho.dims[0]].name: 'time_mid'})
    ds['time_coarse']=(ds.time_coarse/3600)
    ds['time_mid']=(ds.time_mid/3600)
    ds['time_fine']=(ds.time_fine/3600)
    ds['x'] = ds.x.astype(float)*50*1e-3  ### horizontal resolution
    ds['y'] = ds.y.astype(float)*50*1e-3
    return ds


ds_low = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/data/sct_base_files/sct_base_256/sct_base_256_merged.nc')
ds_low = ds_fix_dims(ds_low)

ds_high = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em0_fullA2/sct_em0_merged.nc')
ds_high=ds_fix_dims(ds_high)

fig, ax =plt.subplots(nrows=2, ncols=8, figsize=(20,8))

ds_low.LWP_mean.plot(ax=ax[0,0])
ds_low.RWP_mean.plot(ax=ax[0,0])

ds_low.LWP_mean.plot(ax=ax[0,0])
ds_low.RWP_mean.plot(ax=ax[0,0])


cl_num1=ds_low.q_cloud_liquid_number*1e-6
cl_num1[0,:,128,:].plot(y='z',ax=ax[0,1])

ac_num1=ds_low.q_accum_sol_number*1e-6
ac_num1[0,:,128,:].plot(y='z',ax=ax[0,2])

ait_num1=ds_low.q_aitken_sol_number*1e-6
ait_num1[0,:,128,:].plot(y='z',ax=ax[0,3])

cl_mass1=ds_low.q_cloud_liquid_mass*1e-6
cl_mass1[0,:,128,:].plot(y='z',ax=ax[0,4])

ac_mass1=ds_low.q_accum_sol_mass*1e-6
ac_mass1[0,:,128,:].plot(y='z',ax=ax[0,5])

ait_mass1=ds_low.q_aitken_sol_mass*1e-6
ait_mass1[0,:,128,:].plot(y='z',ax=ax[0,6])

act1=ds_low.q_active_sol_liquid*1e-6
act1[0,:,128,:].plot(y='z',ax=ax[0,7])

cl_num2=ds_high.q_cloud_liquid_number*1e-6
cl_num2[0,:,128,:].plot(y='z',ax=ax[1,1])

ac_num2=ds_high.q_accum_sol_number*1e-6
ac_num2[0,:,128,:].plot(y='z',ax=ax[1,2])

ait_num2=ds_high.q_aitken_sol_number*1e-6
ait_num2[0,:,128,:].plot(y='z',ax=ax[1,3])

cl_mass2=ds_high.q_cloud_liquid_mass*1e-6
cl_mass2[0,:,128,:].plot(y='z',ax=ax[1,4])

ac_mass2=ds_high.q_accum_sol_mass*1e-6
ac_mass2[0,:,128,:].plot(y='z',ax=ax[1,5])

ait_mass2=ds_high.q_aitken_sol_mass*1e-6
ait_mass2[0,:,128,:].plot(y='z',ax=ax[1,6])

act2=ds_high.q_active_sol_liquid*1e-6
act2[0,:,128,:].plot(y='z',ax=ax[1,7])

fig.savefig("aerosol_comp_mass.png")
