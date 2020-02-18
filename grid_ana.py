#!/usr/bin/python

# Description:
# ------------
#     Bla bla bla bla bla
#

import numpy as np
try:
    import matplotlib as mtp
    import matplotlib.pyplot as plt
    from   matplotlib.colors import LinearSegmentedColormap
except ImportError:
    print ("Could not import matplotlib, please make sure you don't")
    print ("want to plot anything.")


# This function plots a generic grid into a matplotlib axis.
#
#   ax <Matplotlib.Axis>:
#       Matplotlib axis instance into which the plot will be carried out
#   title <string>:
#       Title of the grid.
#   x_label <string>:
#       Label for the x axis.
#   y_label <string>:
#       Label for the y axis.
#   cbox_label <string>:
#       Label for the color box.
#   cbox_limits <tuple>:
#       A tuple of floats (m,M) specifying the cut-off values for the color box.
#   cbox_limiting_colors <tuple>:
#       A tuple of strings ('cm', 'cM') or ('cm', 'cmid', 'cM') specifying the colors used in the limits
#       of the color box.
#   data <numpy bi-dimensional array>:
#       A numpy bi-dimensional array containing the data to be plotted.
#       TODO Specify how the grid will be plotted.
#   
def put_grid_into_axis (ax, title, x_label, y_label, cbox_label, cbox_limits, \
        cbox_limits_colors, data):

    plot_data = np.transpose(data)

    # make list of tuples
    if ( (cbox_limits != ()) and (cbox_limits_colors != ()) ):
        if (len(cbox_limits_colors) == 3):
            aux_list = [ (0,cbox_limits_colors[0]), (0.5, cbox_limits_colors[1]),\
                         (1.0, cbox_limits_colors[2]) ]
        if (len(cbox_limits_colors) == 2):
            aux_list = [ (0,cbox_limits_colors[0]), \
                         (1.0, cbox_limits_colors[1]) ]

        # make colormap
        cm = LinearSegmentedColormap.from_list (name='my_cmap', colors=aux_list)
        # plot with colormap
        im = ax.matshow (plot_data, origin='lower', cmap=cm,\
                vmin=cbox_limits[0], vmax=cbox_limits[1])
    else:
        im = ax.matshow (plot_data, origin='lower')

    ax.set_title (title)
    ax.set_xlabel (x_label)
    ax.set_ylabel (y_label)
    ax.xaxis.set_tick_params(width=0)
    ax.set_xticklabels(labels=[])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_tick_params(width=0)
    ax.set_yticklabels(labels=[])
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbox_label, labelpad=18)

    return

# shape = (nx, ny)
def put_markers_into_axis (ax, marker_ids, shape):
    x = []
    y = []
    for m in marker_ids:
        x.append ( int(m) / int(shape[1]) )
        y.append ( int(m) % int(shape[1]) )
    ax.scatter( x,y, color='black' )

# This function plots a grid to a file.
# See function 'put_grid_into_axis' for an overview of the parameters.
def plot_grid_to_file (filename, title, x_label, y_label, cbox_label, cbox_limits, \
        cbox_limits_colors, data, markers):

    mtp.rcParams['text.usetex'] = True
    mtp.rcParams['font.size'] = 20
    mtp.rcParams['font.family'] = 'serif'
    mtp.rcParams['font.serif'] = 'cm'

    fig = plt.figure()
    ax = fig.gca()

    put_grid_into_axis (ax, title, x_label, y_label, cbox_label, cbox_limits, \
        cbox_limits_colors, data)

    put_markers_into_axis (ax, markers, data.shape)

    for axs in fig.get_axes():
        axs.xaxis.set_tick_params(direction='out')
        axs.yaxis.set_tick_params(direction='out')

    plt.tight_layout(True)
    fig.savefig(filename)



# This is only for testing purposes.
if __name__ == '__main__':

    data = np.loadtxt('tmp_data.dat')

    plot_grid_to_file ("grid-image.pdf", "$\\rho_\\mathrm{liq}$", "$x$", "$y$",\
            "$\Delta \\rho_\\mathrm{liq}$", (0,8), ('red','blue') , data)
#def put_grid_into_axis (ax, title, x_label, y_label, cbox_label, cbox_limits, \
#       cbox_limits_colors, x_limits, y_limits, data):

