import re
import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

def save(filename, contents):
    fh = open(filename, 'w')
    fh.write(contents)
    fh.close()


pth = r'F:\Data\Electrodes\MircoExpressions\Experiment\SJ05_A07_B01_ME_stimuli_sj1_block1_6_5_2018_5-42-50PM\Events\sj_6_me.txt'
pat_pth = r'F:\Data\Electrodes\MircoExpressions\Experiment\SJ05_A07_B01_ME_stimuli_sj1_block1_6_5_2018_5-42-50PM\ME.'

file = open(pth, 'r')
strings = file.read()
lines = re.split(r'Event number', strings)
# extracting column headers from the first line
headers = re.split('chan', lines[0]) # keeping only numbers after the file description
lines.pop(0)
# first header is always time
headers.pop(0), headers.pop(0)
headers.insert(0, '\ttime\t')

# adding the relevant headers to the text files
# loading the txt to pandas
for chunk in lines:
    ME_num = chunk.split('\n')[0]
    n_header = (ME_num+'_ch').join(headers)
    f_ind = chunk.find('\t')
    contents = n_header + chunk[f_ind:-1]
    save(pat_pth+ME_num+r'.txt', contents)
# loading txt files to dataFrames
pathlist = Path(pat_pth).glob('**/*.txt')

fig_list = list()
for pt in pathlist:
    df = pd.read_table(str(pt), sep='\t')
    df = df.drop(df.columns[df.columns.str.contains('unnamed', case=False)], axis=1)
    df['time'] = df['time']/1000
    df.plot(x='time', subplots=True, figsize=(15, 15), legend=False, ylim=(-500, 500), color='k')
    ME_num = df.columns.values[1].split('_')[0]
    plt.savefig(pat_pth+ME_num+'.svg', format='svg')

