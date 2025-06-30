import xarray as xr
import sys
sys.path.append("./py_scripts")
import cloud_lib as cl
import re

def create_dataset(member_key, ice=False, count=None):
    sct_dir = "/gws/nopw/j04/carisma/eers/sct"

    if ice:
        splits = member_key.split("_")

        print(f"Opening dataset {member_key}")
        ds = xr.open_dataset(f"{sct_dir}/ice/{splits[0]}/{member_key}/sct_{member_key}_merged.nc")
    else:
        r_key = re.compile(r'(\D*)(\d*)')
        a = r_key.search(member_key)
        key = a.group(1)
    
        print(f"Opening dataset {member_key}")
        ds = xr.open_dataset(f"{sct_dir}/{key}/{member_key}/sct_{member_key}_merged.nc")
    ds = cl.ds_fix_dims(ds)
    
    print('Calculating LWP mask')
    total_lwp, cloudy_lwp, cloud_frac, times, lwp_mask_2d = cl.lwp_cloud(ds)
    
    print('Calculating delta theta threshold')
    deltheta_thresh_heights = cl.deltheta_thresh_timeseries(ds)

    cloud_mass_lim = 1e-5
    rain_mass_lim = 1e-6
    
    print("Calculating in-cloud mean masses")
    cloud_droplet_mass = ds.q_cloud_liquid_mass.where(ds.q_cloud_liquid_mass>cloud_mass_lim).mean(axis=(1,2,3))
    rain_droplet_mass = ds.q_rain_mass.where(ds.q_rain_mass>rain_mass_lim).mean(axis=(1,2,3))
    activated_aerosol_mass = ds.q_cloud_liquid_mass.where(ds.q_active_sol_liquid>cloud_mass_lim).mean(axis=(1,2,3))

    print("Calculating out-of-cloud mean masses")
    aitken_mass = ds.q_aitken_sol_mass.where(ds.q_cloud_liquid_mass<cloud_mass_lim).mean(axis=(1,2,3))
    accumulation_mass = ds.q_accum_sol_mass.where(ds.q_cloud_liquid_mass<cloud_mass_lim).mean(axis=(1,2,3))
    coarse_mass = ds.q_coarse_sol_mass.where(ds.q_cloud_liquid_mass<cloud_mass_lim).mean(axis=(1,2,3))
    
    print("Calculating in-cloud mean numbers")
    cloud_droplet_number = ds.q_cloud_liquid_number.where(ds.q_cloud_liquid_mass>cloud_mass_lim).mean(axis=(1,2,3))
    rain_droplet_number = ds.q_rain_number.where(ds.q_cloud_liquid_mass>cloud_mass_lim).mean(axis=(1,2,3))

    print("Calculating out-of-cloud mean numbers")
    aitken_number = ds.q_aitken_sol_number.where(ds.q_cloud_liquid_mass<cloud_mass_lim).mean(axis=(1,2,3))
    accumulation_number = ds.q_accum_sol_number.where(ds.q_cloud_liquid_mass<cloud_mass_lim).mean(axis=(1,2,3))
    coarse_number = ds.q_coarse_sol_number.where(ds.q_cloud_liquid_mass<cloud_mass_lim).mean(axis=(1,2,3))

    print("Calculating parameter timeseries")
    q_accum_2d_time_z = ds.q_accum_sol_number.mean(axis=(1,2))
    q_cloud_liquid_mass = ds.q_cloud_liquid_mass[:,:,128,:].where(ds.q_cloud_liquid_mass[:,:,128,:]>1e-5)

    rwp = ds.RWP_mean
    toa_down_SW_mean = ds.toa_down_SW_mean
    theta_mean = ds.theta_mean
    vapour_mmr_mean = ds.vapour_mmr_mean
    time_fine = ds.time_fine
    time_mid = ds.time_mid
    time_coarse = ds.time_coarse
    x = ds.x
    y = ds.y
    z = ds.z
    zn = ds.zn
    ds.close()

    print("Creating new xarray dataset")
    post_processed_ds = xr.Dataset({
        "total_lwp":(["time_coarse"], total_lwp),
        "cloudy_lwp":(["time_coarse"], cloudy_lwp),
        "cloud_frac":(["time_coarse"], cloud_frac),
        "lwp_mask_2d":(["time_coarse","x","y"], lwp_mask_2d),
        "inversion_height":(["time_mid"], deltheta_thresh_heights),
        "q_cloud_liquid_mass_3d_masked":(["time_coarse","x","z"], q_cloud_liquid_mass.data),
        "q_accum_sol_number_2d_time_z":(["time_coarse","z"], q_accum_2d_time_z.data),
        "rwp":(["time_fine"], rwp.data),
        "toa_down_SW_mean":(["time_coarse"], toa_down_SW_mean.data),
        "theta_mean":(["time_mid","zn"], theta_mean.data),
        "vapour_mmr_mean":(["time_mid","zn"], vapour_mmr_mean.data),
        "cloud_droplet_mass":(["time_coarse"], cloud_droplet_mass.data),
        "rain_droplet_mass":(["time_coarse"], rain_droplet_mass.data),
        "activated_aerosol_mass":(["time_coarse"], activated_aerosol_mass.data),
        "aitken_mass":(["time_coarse"], aitken_mass.data),
        "accumulation_mass":(["time_coarse"], accumulation_mass.data),
        "coarse_mass":(["time_coarse"], coarse_mass.data),
        "cloud_droplet_number":(["time_coarse"], cloud_droplet_number.data),
        "rain_droplet_number":(["time_coarse"], rain_droplet_number.data),
        "aitken_number":(["time_coarse"], aitken_number.data),
        "accumulation_number":(["time_coarse"], accumulation_number.data),
        "coarse_number":(["time_coarse"], coarse_number.data),
        },
        coords={"time_fine": time_fine, 
                "time_mid": time_mid, 
                "time_coarse": time_coarse,
                "x": x,
                "y": y,
                "z": z,
                "zn": zn
               })

    del lwp_mask_2d, total_lwp, cloudy_lwp, cloud_frac, deltheta_thresh_heights
    for var_ds in [rwp,
                   cloud_droplet_mass,rain_droplet_mass,activated_aerosol_mass,
                   aitken_mass,accumulation_mass,coarse_mass,
                   cloud_droplet_number,rain_droplet_number,
                   aitken_number,accumulation_number,coarse_number,
                   time_fine,time_mid,time_coarse,x,y]:
        var_ds.close()

    print("Saving to netcdf")
    if ice:
        post_processed_ds.to_netcdf(f"{sct_dir}/ice/{splits[0]}/{member_key}/sct_{member_key}_pp.nc",mode="w")
    else:
        post_processed_ds.to_netcdf(f"{sct_dir}/new_processed/sct_em{count}_pp.nc",mode="w")
    post_processed_ds.close()
    print(f"Completed {member_key}!")
    

def main():
    # Process main ensemble
    
    # index_range = [1,4,5,6,8]
    # key_iterator = [f"{key}{index}" for index in index_range]
    # for member_key in key_iterator:
    #     create_dataset(member_key)

    # count = 0
    # key = "em"
    # index_range = range(0,61)
    # key_iterator = [f"{key}{index}" for index in index_range]
    # for member_key in key_iterator:
    #     create_dataset(member_key, count=count)
    #     count += 1

    # key = "val"
    # index_range = range(24)
    # key_iterator = [f"{key}{index}" for index in index_range]
    # for member_key in key_iterator:
    #     create_dataset(member_key, count=count)
    #     count += 1

    # key = "xtra"
    # index_range = range(12)
    # key_iterator = [f"{key}{index}" for index in index_range]
    # for member_key in key_iterator:
    #     create_dataset(member_key, count=count)
    #     count += 1

    key_iterator = [#"em0_0", "em0_1", "em0_2", "em0_3", 
                    #"em34_1", "em34_2", "em34_3", "em34_4",
                    "em80_0", "em80_1", "em80_2", "em80_3", "em80_4",
                   "em85_0", "em85_1", "em85_2", "em85_3", "em85_4"]
    for member_key in key_iterator:
        create_dataset(member_key, ice=True)

if __name__=="__main__":
    main()