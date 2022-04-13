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

font = {'family':'sans-serif', 'size'   : 18}

rc('font', **font)

def ds_fix_dims(ds):
    '''
    Renames time dimensions to fine, mid and coarse and converts times to hours. 
    Converts x and y from number of grid boxes to km. 
    '''
    ### Rename times
    time_dims = [element for element in list(ds.dims.keys()) if "time" in element]
    time_dims.sort(key=lambda s: len(ds[s]))
    ds = ds.rename({ds[time_dims[0]].name:'time_coarse', ds[time_dims[1]].name:'time_mid', ds[time_dims[2]].name: 'time_fine'})
    
    ### Convert to hours
    ds['time_coarse']=(ds.time_coarse/3600)
    ds['time_mid']=(ds.time_mid/3600)
    ds['time_fine']=(ds.time_fine/3600)
    
    ### Convert x and y to km
    dxx_ind = find_options(ds, b'dxx')
    dxx = float(ds.options_database[dxx_ind[0]].values[1])
    ds['x'] = ds.x.astype(float)*dxx*1e-3
    
    dyy_ind = find_options(ds, b'dyy')
    dyy = float(ds.options_database[dyy_ind[0]].values[1])
    ds['y'] = ds.y.astype(float)*dyy*1e-3
    return ds

def smooth_lwp(da, period):
    smoothed_lwp=np.empty(len(da)-period)
    times=np.empty(len(da)-period)
    for i in range(len(da)-period):
        smoothed_lwp[i]=da[i:i+period].mean()
        times[i]=da.time_fine[i:i+period].mean()
    return smoothed_lwp,times

def lwp_cloud(ds):
    lmmr=ds.q_cloud_liquid_mass
    lwp=ds.lwp*1000
    
    total_lwp, cloudy_lwp, cloud_frac, times, lwp_mask = lwp_cloud_calc(lmmr, lwp)
    return total_lwp, cloudy_lwp, cloud_frac, times, lwp_mask

def lwp_cloud_calc(lmmr, lwp):   
    '''
    Calculates the column mask and applies to the lwp with nan where no cloudy gridboxes in the column
    Returns timeseries of total lwp and cloudy lwp, for comparison.
    '''
    total_lwp = np.empty(len(lmmr))
    cloudy_lwp = np.empty(len(lmmr))
    cloud_frac = np.empty(len(lmmr))
    times = np.empty(len(lmmr))
    
    col_mask_array=np.empty(lmmr.shape)
    lwp_mask_array=np.empty(lmmr[:,:,:,0].shape)
    for t in range(len(lmmr)):
        col_mask = layer_cloud_mask_array(lmmr, t, col_mask_array)
        arr_mask = col_mask.mask
        lwp_mask_array[t] = lwp[t].where(arr_mask==True)
        total_lwp[t]=lwp[t].mean(axis=(0,1)).item()
        cloudy_lwp[t]=np.nanmean(lwp_mask_array[t],axis=(0,1))
        times[t]=(lmmr[t][lmmr.dims[0]].item())
        total = arr_mask.sum()
        fraction = total/(len(lmmr.y)*len(lmmr.x))
        cloud_frac[t]=fraction

    return total_lwp, cloudy_lwp, cloud_frac, times, lwp_mask_array


def layer_cloud_mask_array(dataarray, time, col_mask_array):
    '''
    Applies mask to each timestep and sums through the layers
    '''
    for n in range(len(dataarray.z)):
        layer = dataarray[time,:,:,n]
        col_mask_array[time,:,:,n] = np.ma.masked_less(layer, 1e-5, copy=True).mask
    col_sum = np.nansum(col_mask_array[time], axis=2)
    col_mask = np.ma.masked_less(col_sum, len(col_mask_array[time,0,0,:]), copy=True)
    return col_mask


def entrainment(ds):
    '''
    Mean entrainment velocity --- doesn't work
    '''
    ds_qt = ds.vapour_mmr_mean + ds.liquid_mmr_mean + ds.rain_mmr_mean
    h_ind=[]
    [h_ind.append(abs(time - 0.008).argmin().item()) for time in ds_qt]
    #print(ds_qt[np.argmax(h_ind)][np.max(h_ind)].item()-0.008)
    
    divergence = 0.00000375
    we_timeseries = []
    for i,h in enumerate(h_ind[:-1]):
        dz = ds['zn'][h_ind[i+1]].item() - ds['zn'][h].item()
        dt = ds[ds_qt.dims[0]][i+1].item() - ds[ds_qt.dims[0]][i].item()
        we = dz/dt + divergence*ds['zn'][h].item()
        we_timeseries.append(we)
    return we_timeseries

def cloud_bounds_bas_top(ds):
    '''
    calculate cloud boundaries
    '''
    clbas_ave = np.mean(ds.clbas.where(ds['clbas']!=0.0), axis=(1,2))
    cltop_ave = np.mean(ds.cltop.where(ds['cltop']!=0.0), axis=(1,2))
    #print('Summing water...')
#     qt_mass = ds.q_vapour + ds.q_cloud_liquid_mass + ds.q_rain_mass
#     qt_mass_mean = np.mean(qt_mass, axis=(1,2))

#     #t = qt_mass_mean.dims[0]
#     #timeseries=qt_mass_mean[t]
#     #indices = []
    
#     ## find index in each time step closest to isoline
#     #for time in qt_mass_mean:
#     #    idx = np.abs(time - 0.008).argmin()
#     # #   indices.append(int(idx.values))  
#     print('Finding indices...')
#     indices = [int(np.abs(time - 0.008).argmin().values) for time in qt_mass_mean]
        
#     # find corresponding height
#     #height = []
#     #for i in indices:
#     #    height.append(int(qt_mass_mean[qt_mass_mean.dims[1]][i].values))
#     print('Finding heights...')
#     height = [int(qt_mass_mean[qt_mass_mean.dims[1]][i].values) for i in indices]  
    return clbas_ave, cltop_ave #, height

def cloud_bounds_heights(ds):
    '''
    calculate cloud boundaries
    '''
    qt_mass = ds.q_vapour + ds.q_cloud_liquid_mass + ds.q_rain_mass
    qt_mass_mean = np.mean(qt_mass, axis=(1,2))
    print('Finding indices...')
    indices = [int(np.abs(time -0.008).argmin().values) for time in qt_mass_mean]
    print('Finding heights...')
    height=[int(qt_mass_mean[qt_mass_mean.dims[1]][i].values) for i in indices]
    
    return height

# def plot_cb(height, qt_mass_mean, clbas_ave, cltop_ave,fig,ax):
#     '''
#     plot cloud boundaries
#     '''
#     ax.plot(qt_mass_mean[qt_mass_mean.dims[0]], height, color="navy", label="Average inversion height")
#     cltop_ave.plot(color="blue", label="Average cloud top")
#     clbas_ave.plot(color="cyan", label="Average cloud base")
#     ax.set_ylim(0, 1200)
#     ax.set_title("Cloud Boundaries")
#     ax.set_xlabel("Time (hours)")
#     ax.set_ylabel("Height (m)")
#     ax.legend()
#     #plt.show()
#     return ax

def add_diurnal(ds, ax, ylims):
    night_mask = ds.time_coarse.where(ds.toa_down_SW_mean==0.0, drop=True).values
    inds = []
    for n in range(len(night_mask)-1):
        diff = night_mask[n+1]-night_mask[n]
        if diff>8:
            inds.append(n+1)
    splits=np.split(night_mask, inds)
    for s in splits:
        plot_obj = ax.fill_between(s, ylims[0], ylims[1],color='grey', alpha=0.2)
    return plot_obj

def sw_rad(ax, ds, save, xlabel=None, location=None):
    da_list=[ds.flux_up_SW_mean, ds.flux_down_SW_mean, ds.surface_up_SW_mean, ds.surface_down_SW_mean, ds.toa_up_SW_mean, ds.toa_down_SW_mean]
    label_list=["Flux up", "Flux down", "Surface down", "Surface up", "TOA down", "TOA up"]
    ax = plotter(ax, da_list, label_list, ds, ylabel="Domain mean LW", xlabel=xlabel,  ylims=(-1,1400), location=location, diurnal=True)
    if save==True:
        plt.savefig("sw_rad.png")    
    return ax

def lw_rad(ax, ds, save, xlabel=None, location=None):    
    da_list = [ds.flux_up_LW_mean, ds.flux_down_LW_mean, ds.surface_up_LW_mean, ds.surface_down_LW_mean, ds.toa_up_LW_mean]
    label_list = ["Flux up", "Flux down", "Surface down", "Surface up", "TOA up"]
    ax = plotter(ax, da_list, label_list, ds, ylabel="Domain mean LW", xlabel=xlabel,  ylims=(-1,600), location=location, diurnal=True)
    if save==True:
        plt.savefig("lw_rad.png") 
    return ax

def heating_rates(ax, ds, save, lwp_scale, ylims=None, xlabel=None, location=None):
    lwp_scaled = ds.LWP_mean*lwp_scale
    lwp_scaled.plot(ax=ax, label=f"LWP mean (x{lwp_scale})", alpha=0.3)
    da_list = [ds.SW_heating_rate_mean, ds.LW_heating_rate_mean, ds.total_radiative_heating_rate_mean]
    label_list = ["SW", "LW", "Total"]
    ax = plotter(ax, da_list, label_list, ylabel="Heating Rates", xlabel=xlabel, ylims=ylims, location=location, diurnal=True)
    if save==True:
        plt.savefig("heating_rates.png")    
    return ax
    
def plotter(ax, da_list, label_list, ds, xlabel=None, ylabel=None, ylims=None, location=None, diurnal=None):
    for da, label in zip(da_list, label_list):
        da.plot(ax=ax, label=label)    
        
    ax = format_plot(ax, ds, xlabel=xlabel, ylabel=ylabel, ylims=ylims, location=location, diurnal=diurnal)
    return ax

def format_plot(ax, ds, xlabel=None, ylabel=None, ylims=None, location=None, diurnal=None):
    if diurnal==True:
        add_diurnal(ds, ax, ylims)
    if xlabel==True:
        ax.set_xlabel("Time (hours from 8LT)")
    elif xlabel is not None:
        ax.set_xlabel(xlabel)
    else:
        ax.xaxis.set_visible(False)
    ax.set_ylabel(ylabel)
    ax.set_ylim(ylims)
    ax.set_xlim(0,72)
    if location is not None:
        ax.legend(loc=location)
    return ax

def load_lwp_cmap():
    R = 79/255
    G = 79/255
    B = 79/255
    
    cdict = {'red':   [[0.0,  R, R],
                       [0.25,  R*2, R*2],   # the centre colour value is at 0.25 of the normalise vmax
                       [1.0,  1.0, 1.0]],
             'green': [[0.0,  G, G],
                       [0.25,  G*2, G*2],
                       [1.0,  1.0, 1.0]],
             'blue':  [[0.0,  B, B],
                       [0.25,  B*2, B*2],
                       [1.0,  1.0, 1.0]]}
    
    cmap = colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
    norm = colors.Normalize(vmin=0,vmax=700)  # should be same magnitude as lwp
    cmap.set_bad(color=[0/255, 68/255, 94/255])
    
    return cmap, norm

def find_options(ds, b_string):
    '''
    Finds the value of a setting in the options_database. String must be given as bytes type, e.g. b'dxx'
    '''
    results = []
    for i,n in enumerate(ds.options_database.values):
        if b_string in n[0]:
            results.append(i)
    return results