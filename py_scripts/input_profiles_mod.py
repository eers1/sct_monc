#!/usr/bin/env python3.8
import numpy as np

def rand_heights_pert(b_rand_heights, design, row):
    rand_heights = np.asarray(b_rand_heights)
    inversion_height = design[row, 1]
    rand_heights[1] = inversion_height
    rand_heights[2] = inversion_height + 1
    return rand_heights

def heights_profile_pert(full_heights_profile, design, row):
    heights_profile = np.asarray(full_heights_profile)
    inversion = design[row, 1]
    
    ind = np.nonzero(np.where(full_heights_profile - inversion > 0, 1, 0))[0][0]
    heights_profile[ind - 1] = inversion
    heights_profile = np.delete(heights_profile, slice(1, ind-1))
    heights_profile[2] = inversion + 30
    return heights_profile

def theta_profile_pert(b_theta_profile, design, row, heights_profile, th_gradient):
    '''
    Note that the delt is the perturbation but because it's applied over 30m, the difference is a bit bigger
    '''
    delt = design[row, 2]
    inversion = design[row,1]
    
    intercept = (-1)*th_gradient*(b_theta_profile[0] + delt) + inversion
    theta_profile = [(y - intercept)/th_gradient for y in heights_profile]
    theta_profile[:2] = b_theta_profile[:2]
    return np.asarray(theta_profile)

def q_fields_profiles_pert(b_q_fields_profiles, design, row, heights_profile, qv_gradient, length_of_base):
    qv_bl = design[row, 0]*1e-3
    delq = design[row, 3]*1e-3
    na = design[row, 4]
    inversion = design[row,1]
    
    intercept = (-1)*qv_gradient*(qv_bl + delq) + inversion
    qv_profiles = [(y - intercept)/qv_gradient for y in heights_profile]
    qv_profiles[:2] = np.asarray([qv_bl]*2)
    
    accum_bl_num = na*1e6
    accum_bl_mass = find_mass(accum_bl_num)
    
    b_accum_mass = b_q_fields_profiles[3*length_of_base:4*length_of_base]
    b_accum_num = b_q_fields_profiles[4*length_of_base:5*length_of_base]
    
    accum_mass_profile = np.ones(len(heights_profile))*b_accum_mass[-1]
    accum_mass_profile[:2] = np.asarray([accum_bl_mass]*2)
    
    accum_num_profile = np.ones(len(heights_profile))*b_accum_num[-1]
    accum_num_profile[:2] = np.asarray([accum_bl_num]*2)
    
    aitken_mass_profile = b_q_fields_profiles[length_of_base:2*length_of_base]
    aitken_num_profile = b_q_fields_profiles[2*length_of_base:3*length_of_base]
    len_diff = len(aitken_mass_profile)-len(heights_profile)
    if len_diff > 0:
        aitken_mass_profile = np.delete(aitken_mass_profile, slice(((-1)*len_diff)-1,-1))
        aitken_num_profile = np.delete(aitken_num_profile, slice(((-1)*len_diff)-1,-1))
    elif len_diff < 0:
        aitken_mass_profile = np.insert(aitken_mass_profile, 2, np.asarray([aitken_mass_profile[-1]]*abs(len_diff)))
        aitken_num_profile = np.insert(aitken_num_profile, 2, np.asarray([aitken_num_profile[-1]]*abs(len_diff)))
    
    q_fields_profiles = np.hstack((qv_profiles, aitken_mass_profile, aitken_num_profile, accum_mass_profile, accum_num_profile))
    
    return q_fields_profiles, qv_profiles, accum_mass_profile, accum_num_profile

def find_mass(N_mode):
    V = (4*np.pi/3)*(N_mode)*0.5*(0.2e-6)**3*np.exp(9*(np.log(1.5)**2)/2)   #  N_mode should be in # m^-3
    m = V*1500
    return m

def b_aut_pert(design, row, bin_dict):
    b_aut = np.log10(design[row, 5])
    if bin_dict["bin_1"][0] <= b_aut <= bin_dict["bin_1"][1]:
        b_col = bin_dict["bin_1"][2]
    else:
        b_col = bin_dict["bin_2"][2]
    return b_aut, b_col

def base_profiles():
    base_profiles = {}
    base_profiles["brand_heights"] = [0.0, 900.0, 901.0, 3050.01]
    base_profiles["full_heights_profile"] = [10.0, 142.7, 185.7, 228.7, 271.7, 314.7, 357.7, 401.8, 445.9, 490.0, 534.2, 578.3, 623.6, 668.9, 714.2, 759.4, 804.7, 851.7, 898.6, 945.5, 992.5, 1039.4, 1087.6, 1135.9, 1184.1,
                   1232.3, 1280.5, 1330.2, 1380.0, 1429.7, 1479.4, 1529.1, 1580.2, 1631.3, 1682.3, 1733.4, 1784.5, 1837.1, 1889.7, 1942.3, 1994.9, 2047.5, 2101.4, 2155.4, 2209.3, 2263.2, 2317.2, 
                   2372.5, 2427.8, 2483.1, 2538.4, 2593.7, 2651.5, 2709.2, 2766.9, 2824.7, 2882.4, 2940.1, 2997.9, 3055.6, 3100.0]
    base_profiles["heights_profile"] = np.delete(base_profiles["full_heights_profile"], slice(1,18))
    base_profiles["theta_profile"] = [290.969, 290.969, 302.126, 302.394, 302.661, 302.936, 303.211, 303.486, 303.761, 304.036, 304.319, 304.603, 304.886, 305.169, 305.453, 305.744, 306.035, 306.326, 306.617, 306.909, 
                 307.208, 307.508, 307.808, 308.108, 308.408, 308.715, 309.022, 309.33, 309.637, 309.945, 310.26, 310.575, 310.89, 311.206, 311.521, 311.85, 312.179, 312.508, 312.837, 313.167, 
                 313.496, 313.825, 314.154, 314.483]
    base_profiles["q_fields_profiles"] = [0.0104613, 0.0104613, 0.0041191, 0.00410173, 0.00408435, 0.00406649, 0.00404863, 0.00403078, 0.00401292, 0.00399506, 0.00397665, 0.00395824, 0.00393983, 0.00392142, 0.00390302, 
                     0.0038841, 0.00386519, 0.00384628, 0.00382736, 0.00380845, 0.00378897, 0.00376949, 0.00375002, 0.00373054, 0.00371106, 0.00369109, 0.00367112, 0.00365116, 0.00363119, 
                     0.00361122, 0.00359074, 0.00357026, 0.00354978, 0.0035293, 0.00350881, 0.00348744, 0.00346606, 0.00344468, 0.0034233, 0.00340192, 0.00338054, 0.00335916, 0.00333778, 
                     0.00331641, 1.8424807262418176e-11, 1.8424807262418176e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 
                     2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 
                     2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 
                     2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 
                     2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 
                     2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 
                     2.4566409683224232e-11, 2.4566409683224232e-11, 2.4566409683224232e-11, 150000000.0, 150000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 
                     200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 
                     200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 
                     200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 200000000.0, 1.9749950144932982e-09, 
                     1.9749950144932982e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 
                     1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 
                     1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 
                     1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 
                     1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 
                     1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 1.3166633429955323e-09, 
                     1.3166633429955323e-09, 150000000.0, 150000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 
                     100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 
                     100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 
                     100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0, 100000000.0]
    base_profiles["b_parameter"] = -1.79
    return base_profiles