#!/usr/bin/python

import numpy as np
import os 

# assumes that the propfiles are standard Gromacs xvg files
def extract_uncorrelated_frames (xtc, tpr, propfiles, oxtc):
    xtc = os.path.abspath(xtc)
    tpr = os.path.abspath(tpr)
    oxtc = os.path.abspath(oxtc)

    path_of_preffix = '/'.join(oxtc.split('/')[0:-1])
    # create path if it does not exist
    os.system("mkdir -p " + path_of_preffix)

    skips = []
    for propfile in propfiles:
        x = np.loadtxt(propfile, comments=['@','#'], usecols=(1,))
        [ac, skip] = truncated_autocorr_with_skip(x)
        skips.append(skip)

    actual_skip = np.max(skips)

    os.system("echo 0 | gmx trjconv -f %s -s %s -skip %d -o %s" % (xtc,tpr,actual_skip,oxtc))

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

#def skip(ac):
#    a = 0.0
#    nsamples = len(ac)
#    for i in range(len(ac)):
#        a += (1.0 - 1.0 * i / nsamples) * ac[i]
#        if ( np.isnan(a) ):
#            a = 0.0
#        
#    return int (1 + 2.0*a)

#if __name__ == "__main__":
#    import matplotlib.pyplot as plt
#    x = np.loadtxt("energy.xvg")
#    ac = truncated_autocorr(x)
#    plt.figure()
#    plt.plot(ac)
#    plt.show()
#    print skip(ac)

