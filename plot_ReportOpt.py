#!/usr/bin/python

import numpy as np
import os.path
import matplotlib.pyplot as plt
from optparse import OptionParser
import brewer2mpl

from STEM_parsers import parse_reportopt

def plot_reportopt(df):
    epsilon = 1e-8
    Dark2 = brewer2mpl.get_map('Dark2', 'Qualitative', number=3)
    cost_yrng = [0, df['cost'].max()]

    idx_iteration = df['it'].values
    if np.diff(report_opt['it'].values).max() < epsilon:
        idx_iteration = df.index

    # plot cost function
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(idx_iteration,
             df['cost'].values,
             'kv-',
             label='cost fct')
    ax.set_ylim(cost_yrng)
    ax.set_xlim([ -0.5, max(idx_iteration) + 0.5])
    ax.set_xlabel('iteration')
    ax.set_ylabel('cost function (mol OCS m$^{-3}$)$^2$')

    plt.plot(idx_iteration,
             df['misfit'].values,
             'o--',
             color=Dark2.mpl_colors[0],
             label='misfit')

    plt.plot(idx_iteration,
             df['bckg'].values,
             '^--',
             color=Dark2.mpl_colors[1],
             label='background')

    #plot a vertical line for each "new x"
    idx_newx = idx_iteration[np.where(df['task'].values == 'NEW_X')[0]]
    plt.vlines(x=idx_newx,
               ymin=cost_yrng[0],
               ymax=cost_yrng[1],
               linestyles='dotted',
               label='new x')

    plt.legend(loc='best')
    return(fig)


if __name__ == "__main__":

    # parse path to Report.opt
    parser = OptionParser(usage='plot_ReportOpt [options] file_name',
    description='plot cost function for specified Report.opt file.')
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error('incorrect number of arguments')

    # draw the plot
    fname = args[0]
    report_opt = parse_reportopt(fname)
    fig = plot_reportopt(report_opt)
    plt.show()
