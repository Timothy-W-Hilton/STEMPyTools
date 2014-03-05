### produce a contour plot of scaling factors from a STEM t_obs_pred.dat file
import matplotlib
matplotlib.use('Agg')

import os
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import argparse

from STEM_parsers import parse_inputdat, parse_tobspred, parse_STEM_coordinates, get_all_tobspred_fnames
import STEM_vis
import na_map

def draw_plot(run_dir,
              input_dir,
              fname,
              iter,
              ax=None,
              t_str=None,
              cb_axis=True,
              v_rng=(0.0, 10.0),
              cmap=cm.get_cmap('Blues'),
              extend='neither'):
    #parse input.dat
    inputdat_fname = os.path.join( run_dir, 'input.dat')
    inputdat = parse_inputdat(inputdat_fname)
    #parse t_obs_pred.dat emissions factors
    tobspred_fname = os.path.join( run_dir, fname)
    print "plotting contours for " + tobspred_fname
    tobspred = parse_tobspred(tobspred_fname)['emi_fac']

    # translate t_obs_pred emi_fac values into 124 x 124 grid
    [lon,lat,topo] = parse_STEM_coordinates(
        os.path.join(input_dir, 'TOPO-124x124.nc'))
    tobspred_gridded = STEM_vis.coords_to_grid(tobspred['x'].values,
                                               tobspred['y'].values,
                                               tobspred['emi_fac'].values)

    # map the emi_fac values
    if t_str is 'default':
        t_str =  "STEM emi_fac; large slab test, weak priors"
        if iter is not None:
            t_str = t_str + "; iteration {}".format(iter)
    if t_str is None:
        t_str = os.path.basename(fname)
    m_emifac = na_map.NAMapFigure(t_str=t_str,
                                  cb_axis=cb_axis,
                                  map_axis=ax)
    m_emifac.add_ocs_contour_plot(lon,
                                  lat,
                                  tobspred_gridded,
                                  vmin=v_rng[0],
                                  vmax=v_rng[1],
                                  extend=extend,
                                  cbar_t_str='emi_fac',
                                  colorbar_args={'format': '%0.2f'},
                                  cmap=cmap)

    return(m_emifac)


if __name__ == "__main__":

    default_input_dir = os.path.join('/mnt',
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
                        default=os.getcwd(),
                        help=('the STEM run directory. Must ' +
                              'contain t_obs_pred.dat and input.dat'))
    parser.add_argument('-i', '--input_dir',
                        nargs='?',
                        type=str,
                        default=default_input_dir,
                        help=('the STEM input directory. Must ' +
                              'contain TOPO-124x124.nc'))
    parser.add_argument('--t_obs_pred_filename',
                        nargs='?',
                        type=str,
                        default='',
                        dest='top_fname',
                        help=('the t_obs_pred.dat file containing the data '
                              'to be plotted.  Default is the highest-numbered'
                              't_obs_pred_123.dat file in the run directory.' ))
    parser.add_argument('-o', '--outfile',
                        nargs='?',
                        type=str,
                        default=os.path.join(os.getenv('HOME'),
                                             'Plots',
                                             'emi_fac_map.pdf'),
                        dest='outfile',
                        help=('The PDF file to draw to.  '
                              'Defaults to $HOME/Plots/emi_fac_map.pdf'))
    parser.add_argument('--iter',
                        nargs='?',
                        type=int,
                        dest='iter',
                        help=('The optimization iteration for the current ' +
                              't_obs_pred.dat file'))
    parser.add_argument('--vmin',
                        nargs='?',
                        type=float,
                        dest='vmin',
                        default=None,
                        help=('Minimum value for emissions factor colorscale. '
                              'Default is 0.0'))
    parser.add_argument('--vmax',
                        nargs='?',
                        type=float,
                        dest='vmax',
                        default=None,
                        help=('Maximum value for emissions factor colorscale. '
                              'Default is 10.0'))

    args = parser.parse_args()

    if args.top_fname is '':
        #default to the latest t_obs_pred.dat file
        all_fnames = get_all_tobspred_fnames(args.run_dir)
        iter_num = len(all_fnames)
        top_fname = all_fnames[-1]
    else:
        iter_num = args.iter
        top_fname = args.top_fname

    m = draw_plot(args.run_dir,
                  args.input_dir,
                  top_fname,
                  iter=iter_num,
                  v_rng=np.array((args.vmin, args.vmax)),
                  extend='both')

    print 'writing ' + args.outfile
    m.fig.savefig(args.outfile)
