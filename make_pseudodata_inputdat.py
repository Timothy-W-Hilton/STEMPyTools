"""
short script to merge STEM model [OCS] concentrations from a STEM
t_obs_pred.dat file into a STEM input.dat file.  This is useful for
creating psuedo-observations for a STEM inverse test run from STEM
forward run output.
"""

import os
import os.path
import pandas as pd
import argparse

import STEM_parsers

def create_new_input_dat(input_dat_fname,
						 t_obs_pred_fname,
						 outfile_fname,
						 factor):

	input_dat = STEM_parsers.parse_inputdat(input_dat_fname)
	n_rows = input_dat.shape[0]

	t_obs_pred = STEM_parsers.parse_tobspred(t_obs_pred_fname)['ocs_conc']

	input_dat.ix[:, 'COS'] = t_obs_pred.ix[:n_rows, 'mod'].values * factor

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
	print 'wrote ' + outfile_fname

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
        description=("""short script to merge STEM model [OCS] concentrations from a STEM t_obs_pred.dat file into a STEM input.dat file.  This is useful for creating psuedo-observations for a STEM inverse test run from STEM forward run output."""))
	parser.add_argument('--input_dat',
						nargs='?',
						type=str,
						dest='input_dat_fname',
						default=os.path.join(os.getcwd(), 'input.dat'),
						const=os.path.join(os.getcwd(), 'input.dat'),
						help=('full path to the input.dat file. '
                              'The default is $PWD/input.dat'))
	parser.add_argument('--t_obs_pred',
                        nargs='?',
                        type=str,
						dest='t_obs_pred_fname',
                        default=os.path.join(os.getcwd(), 't_obs_pred.dat'),
						const=os.path.join(os.getcwd(), 't_obs_pred.dat'),
                        help=('full path to the t_obs_pred.dat file. '
                              'The default is $PWD/t_obs_pred.dat.'))
	parser.add_argument('-o', '--outfile',
                        nargs='?',
                        type=str,
						dest='outfile',
                        default=os.path.join(os.getcwd(), 'new_input.dat'),
                        const=os.path.join(os.getcwd(), 'new_input.dat'),
                        help=('Full path to the new input.dat file to write. '
                              'The default is $PWD/new_input.dat.'))
	parser.add_argument('-f', '--factor',
                        nargs='?',
                        type=float,
                        default=1.0,
                        help=('factor to multiply t_obs_pred OCS ' +
							  'concentrations by.  The default is 1.0.'))
	args = parser.parse_args()


	create_new_input_dat(args.input_dat_fname,
						 args.t_obs_pred_fname,
						 args.outfile,
						 args.factor)
