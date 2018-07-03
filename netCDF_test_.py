# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 15:11:56 2018

@author: Granoviter
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
import scipy as sp
import xrft


def find_th(ds, type='filt'):
    th_per_subj = dict()
    for key in ds.data_vars.keys():
        data = ds.sel(type=type).data_vars[key].sel(time=0).data
        th = [(mue, std) for mue, std in zip(np.mean(data, axis=1), np.std(data, axis=1))]
        th_per_subj[key] = th
    return th_per_subj


def do_th(ds, th_per_subj, type='filt'):

    res = []
    for key in ds.data_vars.keys():
        th = np.array(th_per_subj[key])
        data = ds.sel(type=type).data_vars[key].sel(time=0).data
        res.append([np.array([row > t[0] + 5*t[1] for row, t in zip(data, th)])])
    return np.squeeze(np.stack(res, axis=0))

def get_subjects_coords(ds):
    return [i for i in ds.data_vars.variables.keys()]

def segment_data(segmented :xr.DataArray, filt_size=1000, min_samp_size=300):
    filter = np.ones((1,filt_size))
    all_sj = []
    for num in segmented['subject'].values:
        all_ch = segmented.sel(subject=num).values
        sig_after_conv = [np.convolve(d, filt_size, mode='same') for d in all_ch]
        # labeling the different segmented areas
        labels = [sp.ndimage.label(row) for row in sig_after_conv]
        # finding the segments sizes
        segment_sizes = [sp.ndimage.sum(mask_row, labeled_row[0], index=range(1, labeled_row[1] + 1))
                         for mask_row, labeled_row in zip(sig_after_conv, labels)]
        # keep only large enough segments, as defined by parameter minimal_sample_size
        relevant_segments = [segment_sizes_row >= min_samp_size for segment_sizes_row
                             in segment_sizes]
        # get slice objects
        slices = [sp.ndimage.find_objects(labels_row[0]) for labels_row
                  in labels]
        # clear un-relevant slices
        slices = [np.array(slice_ch)[relevant_segments_ch] for slice_ch, relevant_segments_ch
                  in zip(slices, relevant_segments)]
    return all_sj

def differential_channel(da, ch):
    return da.sel(subject=[6, 10, 3, 5, 11]) - da.sel(electrode=15, subject=[6, 10, 3, 5, 11])

def apply_rms(data, window_size=1000):

    # creating the window
    window = np.ones(shape=int(np.ceil(window_size)))/np.ceil(window_size)
    # power
    sqr_data = np.power(data, 2)
    # return RMS
    return np.sqrt(np.convolve(sqr_data, window, 'same'))

def create_rms_data(data, electrode=11):
    data_list=[]
    a = differential_channel(data, electrode)
    for sj in a.coords['subject'].values:
        data_list.append(apply_rms(a.sel(subject=sj)))
    return np.stack(data_list, axis=1)



if __name__ == '__main__':

    #pth = 'F:\Data\Electrodes\MircoExpressions\ExperimentME_data_filt_comb.nc'
    pth = 'F:\Data\Electrodes\MircoExpressions\ExperimentME_data_filt_comb.nc'
    ds = xr.open_dataset(pth)
    #da = differential_channel(ds, 11)

    time = ds.__xarray_dataarray_variable__.coords['time'].values
    #rms_da = xr.DataArray(res, dims=('time', 'subject'), coords={'time': time})
    #rms_da.plot.line(row='subject')
    no_rms = differential_channel(ds.__xarray_dataarray_variable__, 11)
    no_rms.plot.line(row='subject')
    no_rms.plot.show()
    plt.show()
    #rms_da.to_netcdf(r'F:\Data\Electrodes\MircoExpressions\rms_1000.nc')
