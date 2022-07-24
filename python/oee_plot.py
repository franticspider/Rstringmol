"""oee_plot.py
   @author: Susan Stepney
   plotting functions, for 
   "On the Open-Endendness of Detecting Open-Endedness"
   by Stepney and Hickinbotham, Artificial Life, 2023"""

import oee_utils
import numpy as np
import matplotlib.pyplot as plt


FIGSIZE = (5, 3)  # uniform size/style for all figs
SCATTERFIGSIZE = (10, 5)  # uniform size/style for all scatter plots
BOXFIGSIZE = (10, 2.5)  # uniform size/style for all box plots
GRIDFIGSIZE = (10, 3)  # uniform size/style for all grids of 8 T=101 plots
STACKEDFIGSIZE = (10, 3)  # uniform size/style for all stacked figs, reduced to half page width
GRIDFIGSIZE12 = (10, 1.5)  # uniform size/style for all grids of 12 narrow plots
plt.rcParams['font.size'] = '12'  # to make labels visible
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.linewidth'] = 0.5


def plot_lines(ys, fname, color='C0', ticksize=None, ylog=False, ymin=False, ymax=False, yticks=None, title=None):
    fig, ax = plt.subplots(1, 1, figsize=FIGSIZE)

    xs = range(0, len(ys))
    if ys.ndim == 1:  # single line
        ax.plot(xs, ys, color)  # default color = C0 = blueish
    else:
        # multiple lines; ensure enough unique colours
        if ys.shape[1] > 10:        # matplotlib >v2.0 has 10 colours, C0..C9
            colors = plt.cm.tab20(np.linspace(0, 1, 20))  # 20 is enough for these plots
            ax.set_prop_cycle('color', colors)
        ax.plot(xs, ys)

    ax.set_xlim(0, 100)
    plt.margins(x=0)    # no gap between y=0 and x axis
    if ylog:
        ax.set_yscale('log')

    if ymin:
        ax.set_ylim(bottom=ymin)
        # ax.spines.bottom.set_position('zero') # use with more up to date matplotlib version
        ax.spines['bottom'].set_position('zero')  # my matplotlib v3.0.3
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.set_ylim(bottom=0)  # ensure always true zero

    if ymax:
        ax.set_ylim(top=ymax)
    if yticks:
        ax.yaxis.set_major_locator(plt.MaxNLocator(yticks))  # set number of yaxis ticks (incl zero)
    ax.margins(y=0)  # no gap between zero and x axis

    if title:
        plt.title(title, loc='left')

    if ticksize:  # override global param
        ax.tick_params(axis='both', which='major', labelsize=ticksize)

    print(fname)
    plt.savefig(fname, bbox_inches='tight')
    plt.close(fig)


def plot_grid(ys_all, fname, ymax=False, yline=False):

    xmax = len(ys_all[0])  # number of record points in the first chart data
    if xmax == 101:  # runs 1-8
        fig, axes = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=GRIDFIGSIZE)
    else:   # runs 9-20
        fig, axes = plt.subplots(nrows=1, ncols=12, sharex=True, sharey=True, figsize=GRIDFIGSIZE12)
        xmax = 41

    ax_lst = axes.flatten()
    # fig.subplots_adjust(wspace=0, hspace=0)   # no space between subplots

    for i, ys in enumerate(ys_all):
        xs = range(0, len(ys))   # number of time points
        ax = ax_lst[i]

        ax.set_xlim(0, xmax - 1)
        plt.margins(x=0)    # no gap between y=0 and x axis

        ax.set_ylim(bottom=0)  # ensure always true zero
        if ymax:
            ax.set_ylim(top=ymax)
        ax.margins(y=0)  # no gap between zero and x axis

        if ys.ndim == 1:  # single line
            ax.plot(xs, ys, lw=0.5)  # default color = blueish
        else:
            # multiple lines; ensure enough unique colours
            if ys.shape[1] > 10:        # matplotlib >v2.0 has 10 colours, C0..C9
                colors = plt.cm.tab20(np.linspace(0, 1, 20))  # 20 is enough for these plots
                ax.set_prop_cycle('color', colors)
            ax.plot(xs, ys, lw=0.5)

        if yline:
            ax.plot([0, xmax], [yline, yline], ':', lw=0.5, color='gray')  # dotted horiz line at height "yline"

    print(fname)
    plt.savefig(fname, bbox_inches='tight')
    plt.close(fig)


def plot_box(ys, fname, ymax=False, yticks=None, ticksize=None):
    fig, ax = plt.subplots(1, 1, figsize=BOXFIGSIZE)

    ax.violinplot(ys, showmedians=False, showextrema=False)
    ax.boxplot(ys, notch=True, showmeans=True, sym='.', widths=0.2, meanprops={
               "marker": ".", "markerfacecolor": "red", "markeredgecolor": "red"})

    if ymax:
        ax.set_ylim(top=ymax)
    if yticks:
        ax.yaxis.set_major_locator(plt.MaxNLocator(yticks))  # set number of yaxis ticks (incl zero)
    if ticksize:  # override global param, where plots shown full width
        ax.tick_params(axis='both', which='major', labelsize=ticksize)

    print(fname)
    plt.savefig(fname, bbox_inches='tight')
    plt.close(fig)


def plot_density(m_ary, fname, y_max):
    fig, ax = plt.subplots(1, 1, figsize=FIGSIZE)
    ax.set_xlim(0, 100)
    # ax.set_ylim(2, 9)

    step = y_max / 30
    yticks = [0, 10, 20, 30]
    ylabels = [round(t * step) for t in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)

    print(fname)
    plt.imshow(m_ary, cmap='cubehelix_r', extent=(0, 101, 0, 30), origin='lower')
    plt.savefig(fname, bbox_inches='tight')

    # plt.show()
    plt.close(fig)


def plot_scatter(xs, ys, ss, fname, xmax=False, ymax=False):
    fig, ax = plt.subplots(1, 1, figsize=SCATTERFIGSIZE)

    ax.scatter(xs, ys, ss, alpha=0.5)

    if xmax:
        ax.set_xlim(right=xmax)
    if ymax:
        ax.set_ylim(top=ymax)

    print(fname)
    plt.savefig(fname, bbox_inches='tight')
    plt.close(fig)


def plot_stacked(ys, fname, ys_line=None, ylog=False, ymax=False):
    # ys are the several lines to stack
    # ys_line is an option separate single line to draw
    fig, ax = plt.subplots(1, 1, figsize=STACKEDFIGSIZE)

    xs = range(0, ys.shape[1])

    # multiple lines; ensure enough unique colours
    if ys.shape[0] > 10:        # matplotlib >v2.0 has 10 colours, C0..C9
        colors = plt.cm.tab20(np.linspace(0, 1, 20))  # 20 is enough for these plots
        ax.set_prop_cycle('color', colors)

    ax.stackplot(xs, ys, baseline='zero', alpha=0.5)
    if ys_line is not None:  # can't use just "if ys_line", because it is an array
        ax.plot(xs, ys_line, 'k', linewidth=1)

    ax.set_xlim(0, 100)
    if ylog:
        ax.set_yscale('log')
    if ymax:
        ax.set_ylim(top=ymax)

    # override global param as stacked plots are reduced in size
    ax.tick_params(axis='both', which='major', labelsize=plt.rcParams['font.size'] * 2)

    print(fname)
    plt.savefig(fname, bbox_inches='tight')
    plt.close(fig)
