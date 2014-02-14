import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from glob import glob
import pdb

from STEM_parsers import parse_inputdat, parse_tobspred

def draw_boxplots(emifac, ax=None):
    """draws boxplots showing the median, mean, and spread of STEM
    emissions factors vs.  STEM optimization iteration."""

    if ax is None:
        fig, ax = plt.subplots(figsize=(10,6))
    else:
        fig = ax.figure
    if type(emifac) is np.ma:
        emifac = emifac.filled(np.nan)

    means = np.ma.mean(emifac, axis=0)
    bx = ax.boxplot(emifac)
    plt.setp(bx['boxes'], color='black')
    plt.setp(bx['medians'], color='red')
    plt.setp(bx['whiskers'], color='black')
    plt.setp(bx['fliers'], color='black', marker='+')

    h_means = ax.scatter(np.arange(emifac.shape[1])+1, means,
                         marker='x',
                         color='red',
                         label='mean value')
    ax.set_xlabel('STEM optimization iteration')
    ax.set_ylabel('emission scaling factor')

    ax.legend((bx['medians'][0], bx['boxes'][0], h_means),
              ('median', '25 & 75 percentile', 'mean'),
              ncol=2,
              numpoints=1,
              scatterpoints=1,
              bbox_to_anchor=(0., 1.02, 1., .102),
              loc='lower left')

    return(fig, ax)

if __name__ == "__main__":

    run_dir = os.path.join('/home',
                            'thilton',
                            'Stem_emi2_onespecies_big_ocssib',
                            'run.TWH_opt_test_large_slab_weak_prior')

    parser = argparse.ArgumentParser(
        description=("plot boxplots of STEM emifac values " +
                     "vs. optimization iteration for a STEM optimization run"))
    parser.add_argument('-r', '--run_dir',
                        nargs='?',
                        type=str,
                        default=run_dir,
                        help=('the STEM run directory. Must ' +
                              'contain t_obs_pred*.dat and input.dat'))
    args = parser.parse_args()

    fnames = glob(os.path.join(args.run_dir, 't_obs_pred*.dat'))
    if fnames:
        emifac_list = [parse_tobspred(f) for f in fnames]
        emifac = [x['emi_fac']['emi_fac'].values for x in emifac_list]
        emifac = np.transpose(np.array(emifac))
        emifac = np.ma.masked_array(emifac, (np.abs(emifac) - 1.0) < 1e-10)

        fig, ax1 = plt.subplots(figsize=(10,6))

        draw_boxplots(emifac, ax=ax1)
        ax1.set_title('"large slab" test optimization -- "weak" priors, 1.0 removed')

        plt.legend(loc='best')

        plt.show()
