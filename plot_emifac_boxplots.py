import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from glob import glob

from STEM_parsers import parse_inputdat, parse_tobspred

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

        means = np.mean(emifac, axis=0)

        fig, ax1 = plt.subplots(figsize=(10,6))
        bx = plt.boxplot(emifac)
        plt.scatter(np.arange(emifac.shape[1])+1, means,
                    marker='*',
                    label='mean value')
        ax1.set_title('"large slab" test optimization -- "weak" priors')
        ax1.set_xlabel('STEM optimization iteration')
        ax1.set_ylabel('emission scaling factor')

        plt.legend(loc='best')
