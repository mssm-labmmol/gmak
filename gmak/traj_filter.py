import gmak.runcmd as runcmd
import numpy as np
import os
from pymbar.timeseries import statisticalInefficiency
from gmak.config import ConfigVariables


def filter_data(data):
    skip = int(statisticalInefficiency(data))
    return skip, data[::skip]


def filter_datafile(datafile, filteredfile, filter_function=None):
    if filter_function is None:
        filter_function = filter_data
    data = np.loadtxt(datafile, comments=["@", "#"], usecols=(-1,))
    skip, fdata = filter_function(data)
    np.savetxt(filteredfile, fdata)
    return skip

