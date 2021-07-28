#!/usr/bin/env python3 

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import warnings

def step2best(stepdir):
    """Receives a directory and returns the indexes of all the points that are
    contain the best within the .95 confidence interval."""
    fn = os.path.join(stepdir, "optimizer_data.dat.cis")
    data = np.loadtxt(fn, comments=("#",))
    best_central = data[0,1]
    similar  = []
    for i, imin, imax in zip(data[:,0], data[:,-2], data[:,-1]):
        if (best_central <= imax) and (best_central >= imin):
            similar.append(int(i))
        else:
            break
    return similar

def step2prop(stepdir):
    """Receives a directory and returns individual matrixes with the
    properties, errors and scores."""
    fn = os.path.join(stepdir, "full_data.dat")
    data = np.loadtxt(fn, comments=("#",))
    avgs = data[:,1:-1:2]
    errs = data[:,2::2]
    scores = data[:,-1]
    return avgs, errs, scores

def grid2extent(griddir, name='parameters_main'):
    """Receives a grid directory and returns the parameter bounding box
    (left, right, bottom, top)."""
    fn = os.path.join(griddir, name + '.dat')
    data = np.loadtxt(fn)
    return (data[0,0], data[-1,0], data[0,1], data[-1,1])

def idxs2par(idxs, griddir, name='parameters_main'):
    """Receives a list of indexes, the grid directory and returns the list of
    parameters corresponding to these indexes."""
    fn = os.path.join(griddir, name + '.dat')
    data = np.loadtxt(fn)
    pars = [data[i,:] for i in idxs]
    return np.array(pars)

def grid2step(griddir):
    """Receives a directory and returns a list of `step` directories within
    it."""
    stepdirs = []
    step = 0
    while(os.path.exists(f"{griddir}/step_{step}")):
        stepdirs.append(f"{griddir}/step_{step}")
        step += 1
    return stepdirs

def base2grid(basedir):
    """Receives a directory and returns a list of `grid' directories within
    it."""
    griddirs = []
    grid = 0
    while(os.path.exists(f"{basedir}/grid_{grid}")):
        griddirs.append(f"{basedir}/grid_{grid}")
        grid += 1
    return griddirs

def plot_grid(fn, data, shape, extent=None, marks=None):
    """Receives a filename, a data matrix and a shape tuple. Reshapes the data
    FOR EACH COLUMN into the given shape and plots the resulting things into
    the given file. If marks is not None, mark the points with the indexes in
    marks. If extent is not None, use it as the bounding box."""
    if data.ndim != 2:
        # warnings.warn("Only 2D plots for now!")
        data = data.reshape((-1, 1))

    fig, axs = plt.subplots(nrows=1, ncols=data.shape[1])

    try:
        for i, ax in enumerate(axs):
            rd = data[:,i].reshape(shape)
            pl = ax.imshow(np.transpose(rd), origin='lower', extent=extent,
                           aspect='auto')
            fig.colorbar(pl, ax=ax)
    except TypeError:
        i = 0
        ax = axs
        rd = data[:,i].reshape(shape)
        pl = ax.imshow(np.transpose(rd), origin='lower', extent=extent,
                       aspect='auto')
        fig.colorbar(pl, ax=ax)

    if marks is not None:
        X = marks[:,0]
        Y = marks[:,1]
        ax.scatter(X, Y, marker='x', color='orange')
        ax.scatter([X[0],],[Y[0],], marker='o', color='red')

    fig.savefig(fn)

if __name__ == "__main__":
    basedir = os.path.abspath(sys.argv[1])
    shape = [int(s) for s in sys.argv[2:]]

    finalsteps = []

    for griddir in base2grid(basedir):
        ext = grid2extent(griddir)
        for stepdir in grid2step(griddir):
            avgs, errs, scores = step2prop(stepdir)
            similar = step2best(stepdir)
            similar = idxs2par(similar, griddir)
            plot_grid(f"{stepdir}/Avgs.pdf", avgs, shape, extent=ext)
            plot_grid(f"{stepdir}/Errs.pdf", errs, shape, extent=ext)
            plot_grid(f"{stepdir}/Scores.pdf", scores, shape, marks=similar,
                      extent=ext)
        finalsteps.append(f"{stepdir}/Scores.pdf")

    print("\n\nmontage -geometry +0+0 %s Scores.pdf" % ' '.join(finalsteps))

