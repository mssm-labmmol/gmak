import runcmd
import numpy as np
import os 
from pymbar.timeseries import statisticalInefficiency
from config import ConfigVariables

# Wrapper for statisticalIneffiency calculation.
# pymbar.timeseries.statisticalInefficiency fails if auto-covariance is zero.
# However, there are cases when the auto-covariance can be zero, e.g. if the series is constant.
# This wrapper prevents errors in these cases.
def wrapperStatisticalInefficiency(data):
    return statisticalInefficiency(data)

# assumes that the propfiles are standard Gromacs xvg files
def extract_uncorrelated_frames (xtc, tpr, propfiles, oxtc, opropfiles, methods):
    # If method = 'mbar', the filtered trajectory and the filtered properties have the same -skip.
    # If method = 'linear', 'nearest' or 'cubic', the filtered trajectory will have the same -skip as in 'mbar', but each property trajectory will be filtered according to its own statistical inefficiency.
    xtc = os.path.abspath(xtc)
    tpr = os.path.abspath(tpr)
    oxtc = os.path.abspath(oxtc)

    path_of_preffix = '/'.join(oxtc.split('/')[0:-1])
    # create path if it does not exist
    runcmd.run("mkdir -p " + path_of_preffix)

    # do nothing if no properties are specified
    if len(propfiles) == 0:
        return

    skips = []
    for propfile in propfiles:
        x = np.loadtxt(propfile, comments=['@','#'], usecols=(1,))
        skip = wrapperStatisticalInefficiency(x)
        skips.append(int(np.ceil(skip)))

    actual_skip = int(np.max(skips))

    runcmd.run("echo 0 | %s trjconv -f %s -s %s -skip %d -o %s" % (ConfigVariables.gmx, xtc,tpr,actual_skip,oxtc))
    
    # also filter the selected properties
    for i, propfile in enumerate(propfiles):
        x = np.loadtxt(propfile, comments=['@','#'], usecols=(0,1,))
        if (methods[i] == 'mbar'):
            x_skipped = x[::actual_skip]
        elif (methods[i] in ['nearest', 'linear', 'cubic', 'gpr', 'empty']):
            x_skipped = x[::skips[i]]
        else:
            raise ValueError("Unknown filtering method ", methods[i])
        np.savetxt(opropfiles[i], x_skipped)

def truncated_autocorr(x):
    lags = range(len(x))
    corr=[1. if l==0 else np.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    flag = False
    for i in range(len(corr)):
        if corr[i] <= 0:
            corr[i] = 0
            flag = True
        if flag:
            corr[i] = 0
    return np.array(corr)

def truncated_autocorr_with_skip(x):
    lags = range(len(x))
    corr=[1. if l==0 else np.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    flag = False
    skip = 0
    for i in range(len(corr)):
        if (corr[i] <= 0) or (np.isnan(corr[i])):
            corr[i] = 0
            flag = True
        if flag:
            corr[i] = 0
        else:
            skip = skip + 1
    return [np.array(corr),skip]

