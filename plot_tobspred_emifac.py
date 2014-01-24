### produce a contour plot of scaling factors from a STEM t_obs_pred.dat file

import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

from STEM_parsers import parse_inputdat, parse_tobspred
import STEM_vis
import na_map

if __name__ == "__main__":

    run_dir = os.path.join('/home',
                            'thilton',
                            'Stem_emi2_onespecies_big_ocssib',
                            'run.TWH_opt_test_large_slab')
    input_dir = os.path.join('/mnt',
                             'home10',
                             'skulkarni',
                             'StemData21Jul2013',
                             'input')

    parser = argparse.ArgumentParser(
        description=("make a contour plot of scaling factors from a " +
                     "STEM t_obs_pred.dat file"))
    parser.add_argument('-r', '--run_dir',
                        nargs='?',
                        type=str,
                        default=run_dir,
                        help=('the STEM run directory. Must ' +
                              'contain t_obs_pred.dat and input.dat'))
    parser.add_argument('-i', '--input_dir',
                        nargs='?',
                        type=str,
                        default=input_dir,
                        help=('the STEM input directory. Must ' +
                              'contain TOPO-124x124.nc'))
    parser.add_argument('-f', '--filename',
                        nargs='?',
                        type=str,
                        default='t_obs_pred.dat',
                        dest='fname',
                        help=('the t_obs_pred.dat file containing the data ' +
                              'to be plotted'))
    args = parser.parse_args()

    #parse input.dat
    inputdat_fname = os.path.join( args.run_dir, 'input.dat')
    inputdat = parse_inputdat(inputdat_fname)
    #parse t_obs_pred.dat emissions factors
    tobspred_fname = os.path.join( args.run_dir, args.fname)
    tobspred = parse_tobspred(tobspred_fname)['emi_fac']

    # translate t_obs_pred emi_fac values into 124 x 124 grid
    [lon,lat,topo] =STEM_vis.parse_STEM_coordinates(
        os.path.join(args.input_dir, 'TOPO-124x124.nc'))
    tobspred_gridded = STEM_vis.coords_to_grid(tobspred['x'].values,
                                               tobspred['y'].values,
                                               tobspred['emi_fac'].values)

    # map the emi_fac values
    t_str =  "STEM emi_fac - 'large slab' test inverse run"
    m_emifac = na_map.NAMapFigure(t_str=t_str,
                                  cb_axis=True)
    m_emifac.add_ocs_contour_plot(lon,
                                  lat,
                                  tobspred_gridded,
                                  cbar_t_str='emi_fac',
                                  colorbar_args={'format': '%0.2f'})

    plt.show()
