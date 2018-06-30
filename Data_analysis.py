from sklearn.decomposition import FastICA
import load_intan_rhd_format as intan
from scipy.signal import firwin, lfilter, iirnotch
from scipy import ndimage, signal
from aux_functions import *

class DataAnalysis:

    def __init__(self, filename):
        try:
            self.filename = filename
            self.data_dict = intan.read_data(filename)
        except:
            print('Unvalid file \ path for *.rhd data')
        [self.data_vec, self.time_sync, self.trigger_vec tor] = self.sync()
        self.sample_rate = self.data_dict['frequency_parameters']['amplifier_sample_rate']
        self.filtered_data = self.filter_data()
        self.rms_data = self.apply_rms()
        self.minimal_sample_size = 300
        # Braking the data to relevant arrays

    def sync(self):

        digital_channel = self.data_dict['board_dig_in_data']
        trigger_vector = sum(digital_channel[channel]*np.power(2, channel-9) for channel in range(9, 16))
        first_index = np.where(trigger_vector == 127)[0].max()
        last_index = np.where(trigger_vector[first_index:-1] == 27+64)[0].min() + first_index
        # Create trigger vector
        trigger_vector = trigger_vector[first_index:last_index]
        # Create time sync_vector
        time_sync = self.data_dict['t_amplifier'][:][first_index:last_index]
        # Create data vector
        data_vec = self.data_dict['amplifier_data'][:,first_index:last_index]
        return [data_vec ,time_sync-time_sync[0], trigger_vector]
        #save sample rate

    def filter_data(self,num_time_samp=512, low_cut_frq=20, high_cut_frq=350, notch=True):

        nyq = 0.5 * self.sample_rate
        bpf = firwin(num_time_samp, [low_cut_frq, high_cut_frq], nyq=nyq, pass_zero=False, scale=False)
        # change this
        if notch is True:
            b, a = iirnotch(50/nyq, 30)
            filt_data = [lfilter(b, a, self.data_vec[i, :]) for i in range(0, 16)]
            b, a = iirnotch(100 / nyq, 30)
            filt_data = [lfilter(b, a, filt_data[i][:]) for i in range(0, 16)]
        return [lfilter(bpf, 1, filt_data[i][:]) for i in range(0, 16)]

    def apply_rms(self, window_size=40):
        # apply rms on the data
        # window size is in ms - converted to samples
        window_size_samp = window_size*(10**-3)*self.sample_rate
        # creating the window
        window = np.ones(shape=int(np.ceil(window_size_samp)))/np.ceil(window_size_samp)
        # power
        sqr_data = [np.power(self.filtered_data[ch],2) for ch in range(0, len(self.filtered_data))]
        # return RMS
        return [np.sqrt(np.convolve(sqr_data[ch], window, 'same'))
                for ch in range(0, len(self.filtered_data))]

    def apply_ica(self):
        # calculate fast ica on the filtered data, return the transformed data & mixing matrix
        ica = FastICA(max_iter=400)
        X = np.array(self.filtered_data)
        return ica.fit_transform(X.T)

    def apply_threshold(self):
        "A naive signal 1d detection using time series threshold"
        # find average and std for each channel
        thresholds = [(mue, std) for mue, std in zip(np.mean(self.rms_data, axis=1), np.std(self.rms_data, axis=1))]
        # applying threshold - derived from the average and standard deviation of the signal
        masks = [row > th[0] + 4*th[1] for row, th in zip(self.rms_data, thresholds)]
        # labeling the different segmented areas
        labels = [ndimage.label(row) for row in masks]
        # finding the segments sizes
        segment_sizes = [ndimage.sum(mask_row, labeled_row[0], index=range(1, labeled_row[1]+1))
                             for mask_row, labeled_row in zip(masks, labels)]
        # keep only large enough segments, as defined by parameter minimal_sample_size
        relevant_segments = [segment_sizes_row >= self.minimal_sample_size for segment_sizes_row
                             in segment_sizes]
        # get slice objects
        slices = [ndimage.find_objects(labels_row[0]) for labels_row
                  in labels]
        # clear unrelevent slices
        self.filt_slices = [np.array(slice_ch)[relevant_segments_ch] for slice_ch, relevant_segments_ch
                       in zip(slices, relevant_segments)]

        return self.filt_slices

    def find_peaks(self, prominence=1.2, width=300, height=70):
        # this function uses scipy find peaks to detect areas of interest
        self.peaks = [signal.find_peaks(rms_ch, prominence=prominence, width=width, height=height) for rms_ch in self.rms_data]
        #self.peaks = [signal.find_peaks_cwt(rms_ch, widths=[300, 500, 1000], min_snr=1.2) for rms_ch in self.filtered_data]

        return self.peaks

    def correlation_single_sample(self):
        pass

    def cluster(self):
        pass

    def slices_to_csv(self):
        filt_slices = self.filt_slices
        d = {'ME_number':list(), 'Type':list(), 'video_time':list(),
             'length':list(), 'electrodes':list(), 'stimulus_number':list(), 'stimulus_time':list()}
        # initiating lists

        video_time = list()
        length = list()
        electrode_num = list()
        stim_num = list()
        cnt = 0
        for ch_number ,ch_data in enumerate(filt_slices):
            for element in ch_data:
                cnt =+ 1
                video_time.append(self.time_sync[element[0].start])
                length.append((element[0].stop-element[0].start)*1/self.sample_rate)
                electrode_num.append(ch_number)
                stim_num.append(self.trigger_vector[element[0].start])

        self.put_it_dataFrame(d, video_time, length, electrode_num, stim_num)

    def peaks_to_csv(self):
        peaks = self.peaks
        d = {'ME_number':list(), 'Type':list(), 'video_time':list(),
             'length':list(), 'electrodes':list(), 'stimulus_number':list(), 'stimulus_time':list()}
        # initiating lists

        video_time = list()
        length = list()
        electrode_num = list()
        stim_num = list()
        for ch_number, ch_data in enumerate(peaks):
            video_time.extend(self.time_sync[ch_data[0]])
            length.extend((ch_data[1]['widths'])*1/self.sample_rate)
            stim_num.extend(self.trigger_vector[ch_data[0]])
            electrode_num.extend([ch_number for x in range(0, len(ch_data[0]))])

        self.put_it_dataFrame(d, video_time, length, electrode_num, stim_num)


    def put_it_dataFrame(self, d, video_time, length, electrode_num, stim_num):
        d['video_time'] = video_time
        d['length'] = length
        d['electrodes'] = electrode_num
        d['stimulus_number'] = stim_num
        d['ME_number'] = [x for x in range(0, len(video_time))]
        d['Type'] = ['UC' for x in range(0, len(video_time))]
        d['stimulus_time'] = ['UD' for x in range(0, len(video_time))]
        sj_str = self.filename[self.filename.find('S')+1:self.filename.find("_")]
        self.df = pd.DataFrame(data=d)
        self.df.to_csv(path_or_buf=os.getcwd()+'\data_analysis'+sj_str+'.csv')

    def plot_electrodes_section(self, start, length, type='Filtered', electrodes=[x for x in range(0, 16)],
                                font_size=12, electrode_numering=True,
                                elec_ylim=[500 for x in range(0, 16)], hide_y_labels=True):

        dat_to_plot = self.select_data_type(d_type='Filtered')
        # number of electrodes to plot
        num_rows = len(electrodes)
        # extracting indices
        s_idx, l_idx = get_indices_from_time(start, length, self.time_sync, self.sample_rate)
        # creating subplots
        fig, sp_ax = plt.subplots(nrows=num_rows, ncols=1, sharex=True)
        # common x and y label
        fig.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axes
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.grid(False)
        # add labels
        plt.xlabel('Time [s]', fontsize=font_size)
        plt.ylabel(r'Amplitude [$\mu$V]', fontsize=font_size)
        # setting svg
        plt.axis('off')
        plt.gca().set_position([0, 0, 1, 1])
        # plotting
        [ax.plot(self.time_sync[s_idx:s_idx+l_idx], dat_to_plot[i][s_idx:s_idx+l_idx]
                 , Color='black', linewidth=0.3) for i, ax in zip(electrodes, sp_ax)]
        # setting y limit
        [ax.set_ylim([-axis_lim, axis_lim]) for ax, axis_lim in zip(sp_ax, elec_ylim)]
        # setting x limit
        sp_ax[-1].set_xlim([start, start+length])
        # int x axis
        sp_ax[-1].xaxis.set_major_locator(MaxNLocator(integer=True))
        # adding electrode number
        if electrode_numering is True:
            electrodes[0] = 'Electrode ' + str(electrodes[0])
            [ax.text(x=start+length-0.5, y=150, s=str(elec), fontsize=10) for ax, elec in zip(sp_ax, electrodes)]
        # remove xtics and yticks labels for all except the last plot
        [ax.xaxis.set_ticks_position('none') for ax in sp_ax[0:-1]]
        [ax.yaxis.set_ticks_position('none') for ax in sp_ax[0:-1]]
        if hide_y_labels is True:
            plt.setp([ax.get_yticklabels() for ax in sp_ax[0:-1]], visible=False)
        # reducing ylabel font size
        plt.setp([ax.get_yticklabels() for ax in sp_ax], fontsize=font_size)
        # remove top and right frames
        [ax.spines['top'].set_visible(False) for ax in sp_ax]
        [ax.spines['right'].set_visible(False) for ax in sp_ax]
        [ax.spines['bottom'].set_visible(False) for ax in sp_ax]
        fig.savefig('time_elec_section' + '.svg')
        plt.show()

    def plot_multiple_signals_across_single_elec(self, signal_list, start=0, length=10,
                    hide_y_labels=True, font_size=12, set_y_lim=True, elec_ylim=500, selected_elec=2):

        fig, sp_ax = create_subplots(num_rows=len(signal_list), font_size=font_size,
                        signal_numering=True, start=start, length=length,
                        electrodes=[x for x in range(0, len(signal_list))], hide_y_labels=hide_y_labels,
                        set_y_lim=set_y_lim, elec_ylim=elec_ylim,
                        labels=('Time [s]', r'Amplitude [$\mu$V]'), first_sig_title='ME')
        dat_to_plot = self.select_data_type(d_type='Filtered')
        idx_list = [get_indices_from_time(sig+start, length, self.time_sync, self.sample_rate) for sig in signal_list]
        # creating x vector
        x_vec = np.linspace(start, start + length, idx_list[0][1])
        [ax.plot(x_vec, dat_to_plot[selected_elec][idx[0]:idx[0]+idx[1]]
                 , Color='black', linewidth=0.3) for idx, ax in zip(idx_list, sp_ax)]
        fig.savefig('time_elec_multiple_sig' + '.svg')
        #plt.show()
        plt.show()

    def unify_signals(self):
        pass

    def select_data_type(self, d_type='Filtered'):
        if d_type is 'Filtered':
            dat_to_plot = self.filtered_data
        elif d_type is 'RMS':
            dat_to_plot = self.rms_data()
        elif d_type is 'ICA':
            dat_to_plot = self.apply_ica()
        return dat_to_plot

    def all_signals_to_data_frame(self):
        col = ['ch_'+str(i) for i in range(0, len(self.filtered_data))]
        filt_array = np.array(self.filtered_data).T
        self.df_filt = pd.DataFrame(data=filt_array, index=self.time_sync, columns=col)
        self.df_filt.to_csv('filtered_data.csv')


if __name__ == '__main__':
    pass