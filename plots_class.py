import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import numpy as np
#matplotlib.style.use('ggplot')
from matplotlib.ticker import MaxNLocator

def plot_16_ch(pth, name):

    df = pd.read_table(pth + '\\' + name+'.txt')
    f, ax = plt.subplots(16, sharex='all', sharey='all', figsize=(2.07, 3.75))
    for col in df.columns.values:
        if col.split('_')[0].strip().isdigit() is False and col != 'time':
            del df[col]
    time = df.pop('time')/1000
    for ix, col in enumerate(df.columns.values):
        if col.split('_')[0].strip().isdigit() is True:
            ax[ix].plot(time, df[col])
            # Fine-tune figure; make subplots close to each other and hide x ticks for
            # all but bottom plot.
            f.subplots_adjust(hspace=0.2, left=0.2)
    # setting line width and color
    plt.setp([a.get_lines() for a in f.axes], linewidth=0.3, color='k')
    # setting ticks
    [a.tick_params(left='on', right='on', top='off', bottom='off', direction='in', labelleft=False) for a in f.axes[:-1]]
    f.axes[-1].tick_params(left='on', right='on', top='off', bottom='off', direction='in')
    # setting spines (frame lines)
    [(a.spines['top'].set_visible(False),  a.spines['bottom'].set_visible(False)) for a in f.axes[:-1]]
    f.axes[-1].spines['top'].set_visible(False)
    # removing labels for all axes except last one
    plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible='off')
    # last plot font
    plt.setp(f.axes[-1].get_yticklabels(), visible='on', fontsize=8)
    # x font
    plt.setp(f.axes[-1].get_xticklabels(), visible='on', fontsize=8)
    # setting axes limits
    [a.set_ylim((-250, 250)) for a in f.axes]
    # major and minor ticks differently
    [a.tick_params(which='major', length=0) for a in f.axes]
    [a.tick_params(which='minor', length=0) for a in f.axes]
    f.axes[-1].tick_params(length=2)
    # showing only
    f.axes[-1].xaxis.set_major_locator(MaxNLocator(integer=True))
    # setting figure labels
    f.text(0.5, 0.04, 'Time [s]', ha='center')
    f.text(0.04, 0.5, r'Amplitude [$\mu$V]', va='center', rotation='vertical')


def plot_frequency():
    pass

def plot_2ch():
    pass


if __name__ == '__main__':
    pth = r'F:\Data\Electrodes\MircoExpressions\Experiment\SJ05_A07_B01_ME_stimuli_sj1_block1_6_5_2018_5-42-50PM'
    name = r'ME. 1'
    plot_16_ch(pth, name)
    #bb = matplotlib.transforms.Bbox(np.array([[0, 0], [50, 50]]))
    plt.gca()
    plt.savefig('test'+r'.eps', papertype='a4', transparent=True)