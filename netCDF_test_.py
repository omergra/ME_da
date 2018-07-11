# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 15:11:56 2018

@author: Granoviter
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import scipy as sp
import scipy.signal as signal
import pandas as pd

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
    return da - da.sel(electrode=ch)

def apply_rms(data, window_size=1000):

    # creating the window
    window = np.ones(shape=int(np.ceil(window_size)))/np.ceil(window_size)
    # power
    sqr_data = np.power(data, 2)
    # return RMS
    return np.sqrt(np.convolve(sqr_data, window, 'same'))

def apply_savgol(data):
    savgol_list = []
    for i in data.coords['electrode'].values:
        savgol_list.append(signal.savgol_filter(data.sel(electrode=i).values, window_length=10000,
                                                polyorder=0, mode='nearest'))
    return savgol_list

def create_rms_data(data):
    data_list = []
    for elec in data.coords['electrode'].values:
        data_list.append(apply_rms(data.sel(electrode=elec)))
    return data_list

# test for one subject
def data_analysis(da):
    # create differential data
    ##da_diff_1 = differential_channel(da.sel(subject=5, electrode=slice(8, 15)), 12)
    ##da_diff_2 = differential_channel(da.sel(subject=5, electrode=slice(0, 7)), 5)
    # data 1 diff will be the final DataArray
    ##comb_arr = da_diff_1.combine_first(da_diff_2)
    # removing mean
    ##array_no_mean = (comb_arr.values.T-np.mean(comb_arr.values, axis=1)).T
    # rectifying data
    rectified = abs(da)
    #comb_arr = xr.DataArray(rectified, dims=('electrode', 'time'))
    # rms data
    #rms_da = create_rms_data(comb_arr)
    low_pass = 2/1500
    b2, a2 = sp.signal.butter(3, low_pass, btype='lowpass')
    env = sp.signal.filtfilt(b2, a2, rectified.values)

    #env_mf = match_filter(data=env, width=1000, sig=300)
    #find_peaks_cwt(env)

    return env, rectified

def find_peaks(data, prominence=10, width=(500, 2000), height=5):
    # this function uses scipy find peaks to detect areas of interest

    return [signal.find_peaks(ch, prominence=prominence, width=width, height=5) for ch in data]

def find_peaks_cwt(data):
    return [sp.signal.find_peaks_cwt(ch, min_snr=10) for ch in data]

def unify_peaks(peaks, tol=0.5, fps=3000):
    tol = tol*fps
    flat_list = [item for sublist in peaks for item in sublist[0]]
    peaks_only = [list(row[0]) for row in peaks]
    flat_list.sort()
    clean_list = []
    prominance = []
    width = []
    for idx, item in enumerate(flat_list):
        for other in flat_list[idx+1:-1]:
            if check_in_range(low=item-tol, high=item+tol, number=other):
                flat_list.remove(other)
                idx = index_2d(peaks_only, other)
            else:
                idx = index_2d(peaks_only, other)
                feature = get_peaks_data(peaks[idx[0]][1], idx[1])
                if item not in clean_list:
                    clean_list.append(item)
                    prominance.append(feature[0])
                    width.append(feature[1])
                break
    return clean_list, prominance, width

def index_2d(myList, v):
    for i, x in enumerate(myList):
        if v in x:
            return (i, x.index(v))

def get_peaks_data(peaks_dict, idx):

    prom = peaks_dict['prominences'][idx]
    width = peaks_dict['widths'][idx]
    return prom, width

def match_filter(data, width, sig):

    return [gaussian_filter(ch, sig, width) for ch in data]

def check_in_range(low, high, number):

    return low <= number <= high

def gaussian_filter(data, sig, width):

    window = sp.signal.general_gaussian(width, p=0.5, sig=sig)
    filtered = signal.fftconvolve(window, data)
    filtered = (np.average(data) / np.average(filtered)) * filtered
    filtered = np.roll(filtered, int(-np.floor(width/2)))
    return filtered

def zscore(x, window):
    # rolling window
    x_rolling = x.rolling(time=window)
    # compute average and std
    s = x_rolling.std().fillna(0)
    m = x_rolling.mean().fillna(0)
    # returning the z-score on the window
    return (x-m)/s, m, s


if __name__ == '__main__':

    pth = 'F:\Data\Electrodes\MircoExpressions\ExperimentME_data_filt_comb.nc'
    ds = xr.open_dataset(pth)
    da = ds.__xarray_dataarray_variable__
    x = da.sel(subject=5, electrode=1)
    output = data_analysis(x)
    env = output[0]
    rectified = output[1]
    env = xr.DataArray(env, dims=['time'], coords=da.coords.get('time').values)
    output = zscore(env, window=6000)
    z_score = output[0]
    m = output[1]
    s = output[2]
    #
    th_pos = m+2.8*s
    #
    th_pos.plot()
    plt.show()


