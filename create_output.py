#!/usr/bin/env python3.8

import xarray as xr
import numpy as np
from os import path
import sys
sys.path.insert(0, "/home/users/eers/sct")
import cloud_func_lib as cfl

def load_lh_design(design):
    print("Loading LH design")
    
    if design=="em":
        lh_design = np.loadtxt("/home/users/eers/sct/lh_design/SCT_EmulatorInputsDesign.csv", delimiter=',', skiprows=1)
    elif design=="val":
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
    
def get_value(ds, output_name, design, i):
    '''
    Switch case to retrieve the output of interest. 
    Should have individual elifs for different time outputs. 
    Should create functions for retrieving data and any processing within the switch case. 
    '''
    if output_name=="lwp_mean_final":
        output_value = ds.LWP_mean.values[-4]*1000 # temporarily -4 to avoid corrupted final values
    elif output_name=="cloud_frac_final":
        cloud_frac_data = np.loadtxt(f"/home/users/eers/sct/lwp_mask_csvs/sct_{design}{i}_cloud_frac.csv")
        output_value = cloud_frac_data[-1]
    else:
        print("Please select a valid output of interest")
        
    return output_value

def create_set(design, output_name):
    print(f"Creating {output_name} output for {design} simulations")
    lh_design = load_lh_design(design)
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

def main(design, output_name):
    output_set = create_set(design, output_name)
    file_name = f"/home/users/eers/sct/output_data/sct_{design}_{output_name}.csv"
    
    print(f"Saving output set as {file_name}")
    np.savetxt(file_name, output_set, delimiter=",")
    
    print("Finished calculation")
    
main(sys.argv[1], sys.argv[2])