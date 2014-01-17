"""
short script to merge STEM model [OCS] concentrations from a STEM
t_obs_pred.dat file nto a STEM input.dat file.  This is useful for
creating psuedo-observations for a STEM inverse test run from STEM
forward run output.
"""

import os.path
import pandas as pd

ROOT_DIR = os.path.join('/mnt',
                        'home10',
                        'thilton',
                        'Stem_emi2_onespecies_big_ocssib')

input_dat_fname = os.path.join(ROOT_DIR, 'run.TWH_fwd_small_slab', 'input.dat')
input_dat = pd.read_csv(input_dat_fname,
                        sep=None,
                        header=None,
                        skiprows=4,
                        names=('t', 'x', 'y', 'z', 'COS'))

n_rows = input_dat.shape[0]

t_obs_pred_fname = os.path.join(ROOT_DIR, 'run.TWH_fwd_small_slab', 't_obs_pred.dat')
t_obs_pred = pd.read_csv(t_obs_pred_fname,
                         delim_whitespace=True,
                         header=None,
                         names=('t', 'obs', 'mod'))

input_dat.ix[:, 'COS'] = t_obs_pred.ix[:n_rows, 'mod'].values

outfile_fname = os.path.join(ROOT_DIR,
                             'run.TWH_opt_test_small_slab',
                             'input.dat')
outfile = open(outfile_fname, 'w')
outfile.write('{}\t1\n'.format(n_rows))
outfile.write('##t, x, y, z, COS (ppbv) | NOAA\n')
outfile.write('0.08 1\n')
outfile.write('1\n')
outfile.close()
input_dat.to_csv(outfile_fname,
                 index=False,
                 header=False,
                 sep='\t',
                 mode='a')
