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
import time
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
def data_analysis(da, window, init_th):
    # create differential data
    # rectifying data
    rectified = abs(da)
    low_pass = 2/1500
    b2, a2 = sp.signal.butter(3, low_pass, btype='lowpass')
    env = sp.signal.filtfilt(b2, a2, rectified.values)
    # returning to da
    env = xr.DataArray(env, dims=['time'], coords={'time': da.coords['time'].values})
    env_r = env.rolling(time=window)
    bias = env_r.min().fillna(0)
    env = env - bias
    env[env < init_th] = 0
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

def find_peaks_with_running_zscore(x, m, s, th, window):

    data_arr = x.values
    mean_arr = m.values
    std_arr = s.values
    binary = np.zeros(shape=mean_arr.shape, dtype=np.float16)
    filt = data_arr.copy()
    std = np.std(filt[0:win])
    mean = np.mean(filt[0:win])
    for ind, elem in enumerate(data_arr):
        if elem - mean > th*std:
            binary[ind] = 1
            if window < ind < len(std_arr)-1:
                filt[ind] = filt[ind-1]
        std = np.std(filt[ind-window+1:ind])
        mean = np.mean(filt[ind-window+1:ind])
    return binary*10




if __name__ == '__main__':
    win = 100
    pth = 'F:\Data\Electrodes\MircoExpressions\ExperimentME_data_filt_comb.nc'
    ds = xr.open_dataset(pth)
    da = ds.__xarray_dataarray_variable__
    x = da.sel(subject=11)
    values = x.values.astype(np.float16).T
    df = pd.DataFrame(values)
    df.to_csv('F:\Data\subject_11.csv')
