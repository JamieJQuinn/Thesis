import math
import matplotlib
import matplotlib.pyplot as plt

from plotting_parameters import *

def latexify(fig_width=None, fig_height=None, columns=1, square=False):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}, optional
    square: boolean,optional
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
        
    if square:
        fig_height = fig_width
        
    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height + 
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\\usepackage{gensymb}',
                                      '\\usepackage{amsmath}'],
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

def create_axes(n_columns=1, axis_side='left'):
    latexify(columns=2)
    fig, axis = plt.subplots()
    remove_spines(axis, axis_side)
    return fig, axis

def save_plot(filename):
    plt.savefig(filename, pad_inches=PAD_INCHES, bbox_inches = 'tight')
    plt.show()
def create_axes(n_columns=1, axis_side='left'):
    latexify(columns=2)
    fig, axis = plt.subplots()
    remove_spines(axis, axis_side)
    return fig, axis

def save_plot(filename):
    plt.savefig(filename, pad_inches=PAD_INCHES, bbox_inches = 'tight')
    plt.show()

def create_axes(n_columns=1, axis_side='left'):
    latexify(columns=2)
    fig, axis = plt.subplots()
    remove_spines(axis, axis_side)
    return fig, axis

def save_plot(filename):
    plt.savefig(filename, pad_inches=PAD_INCHES, bbox_inches = 'tight')
    plt.show()
