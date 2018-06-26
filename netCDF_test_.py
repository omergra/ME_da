# -*- coding: utf-8 -*-
"""
Created on Sun Jun 24 15:11:56 2018

@author: Granoviter
"""
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import xarray as xr

if __name__ == '__main__':
    pth = 'D:\ExperimentME_data_set.nc'
    ds = xr.open_dataset(pth)
    tst = ds.sel(type='filt', electrode=11)
    tst.data_vars['06'].isel(time=0).plot()
    plt.show()
