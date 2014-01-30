#!/Users/tim/Library/Enthought/Canopy_64bit/User/bin/python

import numpy as np
import os.path
import matplotlib.pyplot as plt
import argparse
import brewer2mpl
import pdb

from STEM_parsers import parse_reportopt

def plot_reportopt(df,cost_yrng,title=None):
    epsilon = 1e-8
    Dark2 = brewer2mpl.get_map('Dark2', 'Qualitative', number=3)

    if cost_yrng is None:
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
    if title:
        ax.set_title(title)
        
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

    parser = argparse.ArgumentParser(
        description=("Plot the progresion of the cost " +
                     "function during a STEM optimization"))
    parser.add_argument('filename',
                        nargs='?',
                        type=str,
                        default='./Report.opt',
                        help='a STEM Report.opt file')
    parser.add_argument('-y', '--yrange',
                        nargs=2,
                        type=int,
                        dest='yrng',
                        help='min and max vertical axis values')
    parser.add_argument('-T', '--title',
                        nargs=1,
                        type=str,
                        dest='title',
                        help='title to appear above the plot')
    args = parser.parse_args()

    # draw the plot
    report_opt = parse_reportopt(args.filename)
    fig = plot_reportopt(report_opt, args.yrng, args.title[0])
    plt.show()
