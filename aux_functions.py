import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import glob, os
from Data_analysis import DataAnalysis
import xarray as xr
import pandas as pd
import re
def get_indices_from_time(start, length, time_sync, sample_rate):
    s_idx = np.where(time_sync <= start)[0][-1]
    l_idx = int(length * sample_rate)
    return s_idx, l_idx


def create_subplots(num_rows=16, font_size=12, signal_numering=False, start=0, length=10,
                    electrodes=[x for x in range(0,16)], hide_y_labels=True, set_y_lim=False, elec_ylim=300,
                    labels=('Time [s]', r'Amplitude [$\mu$V]'), first_sig_title='Electrode'):
    fig, sp_ax = plt.subplots(nrows=num_rows, ncols=1, sharex=True)
    # common x and y label
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    # add labels
    plt.xlabel(labels[0], fontsize=font_size)
    plt.ylabel(labels[1], fontsize=font_size)
    # setting y limit
    if type(elec_ylim) is int:
        elec_ylim = [elec_ylim for x in range(0, num_rows)]
    if set_y_lim is True:
        [ax.set_ylim([-axis_lim, axis_lim]) for ax, axis_lim in zip(sp_ax, elec_ylim)]
    # setting x limit
    sp_ax[-1].set_xlim([start, start + length])
    # int x axis
    sp_ax[-1].xaxis.set_major_locator(MaxNLocator(integer=True))
    # adding electrode number
    if signal_numering is True:
        electrodes[0] = first_sig_title + ' ' + str(electrodes[0])
        [ax.text(x=start+length-0.5, y=125, s=str(elec), fontsize=10, color='grey') for ax, elec in zip(sp_ax, electrodes)]
    # remove xtics and yticks labels for all except the last plot
    [ax.xaxis.set_ticks_position('none') for ax in sp_ax[0:-1]]
    [ax.yaxis.set_ticks_position('none') for ax in sp_ax[0:-1]]
    if hide_y_labels is True:
        plt.setp([ax.get_yticklabels() for ax in sp_ax[0:-1]], visible=False)
    # reducing ylabel font size
    plt.setp([ax.get_yticklabels() for ax in sp_ax], fontsize=10)
    # remove top and right frames
    [ax.spines['top'].set_visible(False) for ax in sp_ax]
    [ax.spines['right'].set_visible(False) for ax in sp_ax]

    return fig, sp_ax

def csv_for_data_set(mother_folder='F:\Data\Electrodes\MircoExpressions\Experiment'):
    for dirpath, dirnames, filenames in os.walk(mother_folder):
        for filename in [f for f in filenames if f.endswith(".rhd")]:
            da = DataAnalysis(os.path.join(dirpath, filename))
            da.find_peaks()
            da.peaks_to_csv()

def electrode_over_sub(mother_folder='F:\Data\Electrodes\MircoExpressions\Experiment', elec=11):
    fig, sp_ax = plt.subplots(nrows=7, ncols=1, sharex=True)
    cnt = -1
    for dirpath, dirnames, filenames in os.walk(mother_folder):
        for filename in [f for f in filenames if f.endswith(".rhd")]:
            da = DataAnalysis(os.path.join(dirpath, filename))
            sp_ax[cnt].plot(da.time_sync, da.filtered_data[elec])
        cnt = cnt+1
    plt.show()



def create_data_set(mother_folder='F:\Data\Electrodes\MircoExpressions\Experiment'):
    "crates an xarray data frame containing all the subjects data"
    sj_number = []
    elec_num = [i for i in range(16)]  # change for parameter
    data_type = ['raw', 'filt', 'rms']
    all_variables = []
    # creating a data array list
    filt_sj_data = []
    raw_sj_data = []
    rms_sj_data = []
    # loading the data
    for dirpath, dirnames, filenames in os.walk(mother_folder):
        for filename in [f for f in filenames if f.endswith(".rhd")]:
            # subject number
            sj_number.append(filename[filename.find('S') + 1:filename.find("_")])
            # create data analysis class
            da = DataAnalysis(os.path.join(dirpath, filename))
            filt = np.array(da.filtered_data)
            rms = np.array(da.rms_data)
            raw = np.array(da.data_vec)
            filt_sj_data.append(filt)
            raw_sj_data.append(raw)
            rms_sj_data.append(rms)

    # convert to data sets
    filt_sj_data, time_sync = cut_array_to_minimum_shape(filt_sj_data, da.time_sync)
    #raw_sj_data, _ = cut_array_to_minimum_shape(raw_sj_data, da.time_sync)
    #rms_sj_data, time_sync = cut_array_to_minimum_shape(rms_sj_data, da.time_sync)
    # concatenating the data
    #filt_sj_data = [data for data in filt_sj_data]

    #raw_sj_data = [data for data in raw_sj_data]
    #rms_sj_data = [data for data in rms_sj_data]
    filt_da = np.stack(filt_sj_data, axis=2)
    #raw_da = np.stack(raw_sj_data, axis=2)
    #rms_da = np.stack(rms_sj_data, axis=2)
    sj_number = [int(re.findall(r'\d+', sj)[0]) for sj in sj_number]
    # creating data set
    data_xr = xr.DataArray(filt_da, dims=('electrode', 'time', 'subject'), coords={'electrode': elec_num,
                                                                                  'time': time_sync,
                                                                                  'subject': sj_number})
    data_xr.to_netcdf(path=mother_folder+'ME_data_filt_comb.nc')


def cut_array_to_minimum_shape(all_sj_data, time_sync):
    # get the desired array size for entire data set
    array_length = min([data.shape[1] for data in all_sj_data])
    corrected_length_all_sj_data = [data[:, 0:array_length] for data in all_sj_data]
    return corrected_length_all_sj_data, time_sync[0:array_length]


def da_dictionaries(all_sj_data, sj_number):
    # creating dictionary for Data set
    new_dict = dict()
    for num, da_arr in zip(sj_number, all_sj_data):
        new_dict[num] = (['electrode', 'time', 'subject'], da_arr)
    return new_dict

def arrange_data(data, sj_vec):
    org_data = []
    for sj in data:
        org_data.append([row for row in sj])
        col = []
    return org_data

if __name__ == '__main__':
    create_data_set()