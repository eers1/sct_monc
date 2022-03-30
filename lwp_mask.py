#!/usr/bin/env python3.8

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
import matplotlib.animation as animation
import numpy as np
from IPython.display import HTML, display, Video
import matplotlib.patches as patches
from copy import copy
from matplotlib.gridspec import GridSpec
import re
import sys
sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl

font = {'family':'sans-serif', 'size'   : 18}

rc('font', **font)

path = sys.argv[1]

r_path = re.compile(r'carisma/eers/(\S*)/(\S*)/(\D*)(\d*)/(\S*)_merged.nc')
a = r_path.search(path)
project = a.group(1)
design = a.group(2)
key = a.group(3) + a.group(4)
run_name = a.group(5)
print(project, design, key, run_name)

print('Opening dataset')
ds = xr.open_dataset(path)
ds = cfl.ds_fix_dims(ds)
print('Calculating LWP mask')
total_lwp, cloudy_lwp, cloud_frac, times, lwp_mask = cfl.lwp_cloud(ds)
print('Calculating cloud boundaries')
height = cfl.cloud_bounds_heights(ds) 
print('Calculating CDNC')
cdnc = ds.q_cloud_liquid_number
cdnc_where = cdnc.where(cdnc>20)
cdnc_mean = cdnc_where.mean(axis=(1,2,3))*1e-6

print('Saving CSVs')
np.savetxt(f'lwp_mask_csvs/{run_name}_total_lwp.csv', total_lwp, delimiter=',')
np.savetxt(f'lwp_mask_csvs/{run_name}_cloudy_lwp.csv', cloudy_lwp, delimiter=',')
np.savetxt(f'lwp_mask_csvs/{run_name}_cloud_frac.csv', cloud_frac, delimiter=',')
np.savetxt(f'lwp_mask_csvs/{run_name}_times.csv', times, delimiter=',')

arr = lwp_mask
arrReshaped = arr.reshape(arr.shape[0], -1)
np.savetxt(f'lwp_mask_csvs/{run_name}_lwp_mask_2d.csv', arrReshaped)
np.savetxt(f'lwp_mask_csvs/{run_name}_inv_height.csv', height, delimiter=',')
np.savetxt(f'lwp_mask_csvs/{run_name}_active_aero.csv', cdnc_mean, delimiter=',')

print('Finished creating LWP masks')