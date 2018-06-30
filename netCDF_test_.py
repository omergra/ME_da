# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 15:11:56 2018

@author: Granoviter
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
from scipy import ndimage

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
        labels = [ndimage.label(row) for row in sig_after_conv]
        # finding the segments sizes
        segment_sizes = [ndimage.sum(mask_row, labeled_row[0], index=range(1, labeled_row[1] + 1))
                         for mask_row, labeled_row in zip(sig_after_conv, labels)]
        # keep only large enough segments, as defined by parameter minimal_sample_size
        relevant_segments = [segment_sizes_row >= min_samp_size for segment_sizes_row
                             in segment_sizes]
        # get slice objects
        slices = [ndimage.find_objects(labels_row[0]) for labels_row
                  in labels]
        # clear un-relevant slices
        slices = [np.array(slice_ch)[relevant_segments_ch] for slice_ch, relevant_segments_ch
                  in zip(slices, relevant_segments)]
    return all_sj

if __name__ == '__main__':

    pth = 'D:\ExperimentME_data_set.nc'
    ds = xr.open_dataset(pth)
    a = find_th(ds, type='filt')
    data_after_threshold = do_th(ds, a, type='filt')
    segmented = xr.DataArray(data_after_threshold, dims=['subject', 'electrode', 'time'], coords={'subject': [int(i.split('J')[-1]) for i in get_subjects_coords(ds)] ,
                                                                                      'electrode': [i for i in range(16)],
                                                                                     'time': ds.coords['voltage'].values})
    all_sj = segment_data(segmented)
    segmented.to_netcdf()
    #segmented.sel(electrode=slice(9,15)).plot.line(row='subject', col='electrode')
    plt.show()
    tst = ds.sel(type='filt', electrode=slice(1, 2), voltage=slice(100, 200))

