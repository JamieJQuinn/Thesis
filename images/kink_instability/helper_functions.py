from plotting_parameters import *
import math
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sdf

sys.path.insert(0,'../lare3d-tools')
from print_energy import Energy

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2,3])

    if fig_width is None:
        if columns == 1:
            fig_width = COLUMN_WIDTH
        elif columns == 2:
            fig_width = COLUMN_WIDTH * COLUMN_HALFSIZE  
        else:
            fig_width = COLUMN_WIDTH * COLUMN_THIRDSIZE

    if fig_height is None:
        golden_mean = (math.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches
        
    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height + 
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\\usepackage{gensymb}', r"\usepackage{amsmath}"],
              'axes.labelsize': FONTSIZE, # fontsize for x and y labels (was 10)
              'axes.titlesize': FONTSIZE,
              'font.size': FONTSIZE, # was 10
              'legend.fontsize': FONTSIZE, # was 10
              'xtick.labelsize': FONTSIZE,
              'ytick.labelsize': FONTSIZE,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)

def remove_spines(axis, axis_side='left'):
    axis.spines['top'].set_visible(False)
    axis.xaxis.set_ticks_position('bottom')
    if axis_side == 'left':
        axis.spines['right'].set_visible(False)
        axis.yaxis.set_ticks_position('left')
        axis.yaxis.set_label_position('left')
    else:
        axis.spines['left'].set_visible(False)
        axis.yaxis.set_ticks_position('right')
        axis.yaxis.set_label_position('right')

def delta(connectivities, j):
    return ((connectivities[j] - connectivities[j-1]) != 0).astype(int)\
#          + ((connectivities[j] - connectivities[j-1]) < 0).astype(int)

def fix_current_restart(data):
    minimum = data[1:].argmin()
    if data[minimum+1] < 0.001:
        data[minimum+1] = 0.5*(data[minimum] + data[minimum+2])
        
def fix_ohmic_heating_restart(data):
    for i in range(2):
        diff = data[1:] - data[:-1]
        if np.any(np.abs(diff) > 0.01):
            max_point = np.argmax(np.abs(diff))
#             print(max_point, data[max_point-1:max_point+1], diff[max_point-1:max_point+1])
            if diff[max_point] > 0:
                data[max_point+1:] += diff[max_point]
            else:
                data[max_point+1:] -= diff[max_point]

def get_maxes(filenames, index, min_time=-1, max_time=-1):
    energies = [Energy(f) for f in filenames]
    return [get_max(en.data, index, min_time, max_time) for en in energies]

def get_max(data, index, min_time=-1, max_time=-1):
    if max_time > 0:
        max_index = find_nearest(data[:,0], max_time)
    else:
        max_index = len(data[:,0])+1
        
    if min_time > 0:
        min_index = find_nearest(data[:,0], min_time)
    else:
        min_index = 0
    
    return np.max(data[min_index:max_index, index])

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def create_axes(n_columns=1, axis_side='left'):
    latexify(columns=2)
    fig, axis = plt.subplots()
    remove_spines(axis, axis_side)
    return fig, axis

def save_plot(filename):
    plt.savefig(filename, pad_inches=PAD_INCHES, bbox_inches = 'tight')
    plt.show()
