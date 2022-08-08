#!/usr/bin/env python3.8

import xarray as xr
import numpy as np
from os import path
import sys
sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl

def load_lh_design(design, post_spinup):
    '''
    Skiprow in loadtxt is to miss heading row
    '''
    
    print("Loading LH design")
    
    if design=="em":
        if post_spinup=="True":
            lh_design = np.loadtxt("/home/users/eers/sct/lh_design/post_spinupvalues/ppe_post_spinup.csv", delimiter=',')
        else:
            lh_design = np.loadtxt("/home/users/eers/sct/lh_design/SCT_EmulatorInputsDesign.csv", delimiter=',', skiprows=1)
    elif design=="val":
        if post_spinup=="True":
            lh_design = np.loadtxt("/home/users/eers/sct/lh_design/post_spinupvalues/val_post_spinup.csv", delimiter=',')
        else:
            lh_design = np.loadtxt("/home/users/eers/sct/lh_design/SCT_ValidationInputsDesign.csv", delimiter=',', skiprows=1)
    elif design=="base":
        lh_design = np.loadtxt("/home/users/eers/sct/lh_design/SCT_BaseInputs.csv", delimiter=',', skiprows=1)
    elif design=="oat":
        lh_design = None
        print("OAT design not configured yet.")
    else:
        lh_design = None
        print("Please select a valid design")
        
    return lh_design

def open_dataset(path):
    ds = xr.open_dataset(path)
    ds = cfl.ds_fix_dims(ds)
    return ds

def get_sc_index(design, i):
    cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
    times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
    if any(t > 0.9 for t in cf):
        sc, = np.where(cf>0.9)
        sc_index=sc[0]
    else:
        sc_index=-1
    return sc_index
    
def get_value(ds, output_name, design, i):
    '''
    Switch case to retrieve the output of interest. 
    Should have individual elifs for different time outputs. 
    Should create functions for retrieving data and any processing within the switch case. 
    '''
    if output_name=="lwp_mean_final":
        output_value = ds.LWP_mean.values[-4]*1000 # temporarily -4 to avoid corrupted final values
    elif output_name=="lwp_mean_init":
        output_value = ds.LWP_mean.values[5]*1000 # starting at hour 3 to allow for spin up
    elif output_name=="cloud_frac_init":
        cloud_frac_data = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv")
        output_value = cloud_frac_data[1] # starting at hour 3 to allow for spin up
    elif output_name=="cloud_frac_final":
        cloud_frac_data = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv")
        output_value = cloud_frac_data[-1]
    elif output_name=="cloud_frac_day2":
        cloud_frac_data = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv")
        output_value = cloud_frac_data[13] # second day starts at 21 hours
    elif output_name=="inversion_final":
        inv=np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_inv_height.csv")
        output_value = inv[-1]
    elif output_name=="cltop_ave_final":
        output_value = np.mean(ds.cltop.where(ds['cltop']!=0.0), axis=(1,2))[-1]
    elif output_name=="transition_time":
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            cu, = np.where(cf[sc[0]:] < 0.55) # find and save the indices after the CF reaches 0.8 where it then falls below 0.5
            if len(cu)!=0 and all(cf[sc[0]+cu[0]:] < 0.9):
                output_value = times[sc[0]+cu[0]] - times[sc[0]]
            else:
                output_value = 80
        else:
            output_value = -1
#     elif output_name=="sct_time_cut":
#         sct_time = get_value(ds, 'sct_time', design, i)
#         sct_time_cut=sct_time[sct_time<80]
#         output_value=sct_time_cut[sct_time_cut>-1]
    elif output_name=="confirmed_transition_time":
        transition_times = np.loadtxt(f"/home/users/eers/sct/output_data/sct_all_confirmed_transition_times_raw.csv", delimiter=',')
        if design=='em':
            output_value = transition_times[i]
        elif design=='val':
            output_value = transition_times[61+i]
        else:
            print('Must use design = all for confirmed_transition_time')
    elif output_name=="rain_mass":
        output_value = ds.q_rain_mass.sum()
    elif output_name=="rain_number":
        output_value = ds.q_rain_number.sum()
    elif output_name=="surface_precip":
        output_value = ds.surface_precip.sum()
    elif output_name=="inv_rate":
        inv_heights = np.loadtxt(f'/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_inv_height_deltheta.csv', delimiter=',')
        gradient_array = np.gradient(inv_heights)
        output_value = np.mean(gradient_array)
    elif output_name=="inv_rate_sc": # after Sc
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            inv_heights = np.loadtxt(f'/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_inv_height_deltheta.csv', delimiter=',')
            gradient_array = np.gradient(inv_heights[sc[0]:])
            output_value = np.mean(gradient_array)
        else:
            output_value = -1
    elif output_name=="first_rain": # after Sc
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            rain, = np.where(np.asarray([t.sum() for t in ds.q_rain_mass[sc[0]:]])!=0)
            if len(rain)!=0:
                output_value = times[sc[0]+rain[0]] - times[sc[0]]
            else:
                output_value = 80
        else:
            output_value = -1
    elif output_name=="rain25": # after Sc
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            rain, = np.where(np.asarray([t.sum() for t in ds.q_rain_mass[sc[0]:]])>25)
            if len(rain)!=0:
                output_value = times[sc[0]+rain[0]] - times[sc[0]]
            else:
                output_value = 80
        else:
            output_value = -1
    elif output_name=="rain50": # after Sc
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            rain, = np.where(np.asarray([t.sum() for t in ds.q_rain_mass[sc[0]:]])>50)
            if len(rain)!=0:
                output_value = times[sc[0]+rain[0]] - times[sc[0]]
            else:
                output_value = 80
        else:
            output_value = -1
    elif output_name=="acc_minSc":
        '''
        Check for first Sc index and find vertical acceleration minimum after that
        '''
        acc_mean = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_accelerations_min.csv", delimiter=',')
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            output_value = np.min(acc_mean[sc[0]:])
        else:
            output_value = -999
    elif output_name=="acc_minSc_time":
        '''
        Check for first Sc index and find vertical acceleration minimum after that and take the time
        '''
        acc_mean = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_accelerations_min.csv", delimiter=',')
        cf = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv", delimiter=',')
        times = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_times.csv", delimiter=',')
        if any(t > 0.9 for t in cf):
            sc, = np.where(cf > 0.9)  # find and save the indices where the CF is over 0.8
            min_acc = np.argmin(acc_mean[sc[0]:])
            output_value = times[sc[0]+min_acc] - times[sc[0]]
        else:
            output_value = -999
    elif output_name=="acc_min":
        '''
        Find vertical acceleration minimum for the whole timeseries
        '''
        acc_mean = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_accelerations_min.csv", delimiter=',')
        output_value = np.min(acc_mean)
    else:
        print("Please select a valid output of interest")
        
    return output_value

def create_set(design, output_name, post_spinup):
    print(f"Creating {output_name} output for {design} simulations")
    lh_design = load_lh_design(design, post_spinup)
    lh_size = np.shape(lh_design)[0]
    output_set = np.c_[lh_design, np.empty(lh_size)]
    
    for i in range(lh_size):
        print(f"Calculating simulation {design}{i}")
        
        ds_path = f"/gws/nopw/j04/carisma/eers/sct/{design}/{design}{i}/sct_{design}{i}_merged.nc"
        file_bool = path.isfile(ds_path)
        
        if file_bool==True:
            ds = open_dataset(ds_path)
            output_set[i,-1] = get_value(ds, output_name, design, i)
            ds.close()
        else:
            print(f"No merged file available, setting {design}{i} output to NaN")
            output_set[i,-1] = np.nan
        
    return output_set

def main(design, output_name, post_spinup):
    if design=="all":
        em_output_set = create_set("em", output_name, post_spinup)
        val_output_set = create_set("val", output_name, post_spinup)
        output_set = np.concatenate((em_output_set, val_output_set), axis=0)
    else:
        output_set = create_set(design, output_name, post_spinup)
    file_name = f"/home/users/eers/sct/output_data/sct_{design}_{output_name}_post_spin_{post_spinup}.csv"
    #file_name = f"/home/users/eers/sct/output_data/sct_{design}_{output_name}.csv"
    
    print(f"Saving output set as {file_name}")
    np.savetxt(file_name, output_set, delimiter=",")
    
    print("Finished calculation")
    
main(sys.argv[1], sys.argv[2], sys.argv[3])