#!/usr/bin/env python3.8

import xarray as xr
import matplotlib.pyplot as plt


def ds_fix_dims(ds):
    ds = ds.rename({ds[ds.w.dims[0]].name:'time_coarse', ds[ds.w_max.dims[0]].name:'time_fine',
	ds[ds.rho.dims[0]].name: 'time_mid'})
    ds['time_coarse']=(ds.time_coarse/3600)
    ds['time_mid']=(ds.time_mid/3600)
    ds['time_fine']=(ds.time_fine/3600)
    ds['x'] = ds.x.astype(float)*50*1e-3  ### horizontal resolution
    ds['y'] = ds.y.astype(float)*50*1e-3
    return ds

DS = {}
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/test_runs/lhb_tests/lhb_micro_3day/sct_lhb_micro_merged.nc')
DS['sct_lhb_micro']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em0_fullA2/sct_em0_merged.nc')
DS['base']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em1_fullA2/sct_em1_merged.nc')
DS['em1']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em2/sct_em2_merged.nc')
DS['em2']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em3/sct_em3_merged.nc')
DS['em3']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em4/sct_em4_merged.nc')
DS['em4']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em5/sct_em5_merged.nc')
DS['em5']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em6/sct_em6_merged.nc')
DS['em6']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em7/sct_em7_merged.nc')
DS['em7']=ds_fix_dims(ds)
ds = xr.open_dataset('/gws/nopw/j04/carisma/eers/sct_sim/ppe/em8/sct_em8_merged.nc')
DS['em8']=ds_fix_dims(ds)

fig,ax = plt.subplots()
for key,val in DS.items():#
    da=val.q_cloud_liquid_number
    da=da.where(da>0)
    da_mean=da.mean(axis=(1,2,3))
    da_mean.plot(ax=ax, label=key)
    
plt.legend()
plt.savefig("cloud_drop_no.png")
