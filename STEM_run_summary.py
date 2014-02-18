import matplotlib
matplotlib.use('Agg')   # this seems to get around the 'couldn't connect to display "localhost:10.0" error' that crops up now and then -- TWH  

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
import os
import os.path
import pdb
import sys
from datetime import datetime
import argparse

import STEM_parsers
import STEM_vis
from na_map import NAMapFigure
import plot_ReportOpt
import plot_emifac_boxplots
import plot_tobspred_emifac

def initialize_plotting_objects(n_plots=8, figsize=(8.5, 11)):
    """create a figure and axes instances to hold STEM run summary
    plots.  Returns a tuple containing the figure and a list of axes"""

    n_rows = np.int(np.ceil(n_plots/2))
    n_cols = 10

    fig = plt.figure(figsize=figsize)
    #define spacing for left-hand column of panels

    gs0 = gridspec.GridSpec(n_rows, n_cols)
    ax_list = {'text_panel':fig.add_subplot(gs0[0, :])}

    gs1 = gridspec.GridSpec(n_rows, n_cols)
    gs1.update(left=0.10, right=0.50, wspace=0.05, hspace=0.30)
    ax_list.update({'cost_func':fig.add_subplot(gs1[1, 0:9]),
                    'emi_fac_map_final':fig.add_subplot(gs1[2, 0:7]),
                    'emi_fac_map_final_cbar':fig.add_subplot(gs1[2, 8]),
                    'pseudo_map':fig.add_subplot(gs1[3, 0:7]),
                    'pseudo_map_cbar':fig.add_subplot(gs1[3, 8])})
    #define spacing for right-hand column of panels
    gs2 = gridspec.GridSpec(n_rows, n_cols)
    gs2.update(left=0.55, right=0.95, wspace=0.05, hspace=0.30)
    ax_list.update({'emi_fac_boxplots':fig.add_subplot(gs2[1, 0:9]),
                    'emi_fac_map_N':fig.add_subplot(gs2[2, 0:7]),
                    'emi_fac_map_N_cbar':fig.add_subplot(gs2[2, 8]),
                    'fwd_map':fig.add_subplot(gs2[3, 0:7]),
                    'fwd_map_cbar':fig.add_subplot(gs2[3, 8])})

    return((fig, ax_list))

def summary_text_panel(ax, run_dir, input_dir, fontsize=8, note=''):
    ax.clear()
    txt = ('STEM run summary - produced ' +
           datetime.strftime(datetime.now(), '%d %b %Y %H:%M:%S') + '\n' +
           'run directory: ' + run_dir + '\n'
           'input directory: ' + input_dir + '\n' +
           note + '\n')
    text_obj = ax.text(0.0,
                       1.0,
                       txt,
                       fontsize=fontsize,
                       horizontalalignment='left',
                       verticalalignment='top',
                       transform = ax.transAxes)
    ax.set_axis_off()
    return(text_obj)

def plot_inputdat_conc(input_dir,
                       run_dir,
                       map_axis,
                       cb_axis,
                       t_str='',
                       cbar_t_str=''):

    input_dat = STEM_parsers.parse_inputdat(os.path.join(run_dir, 'input.dat'))
    gridded_input_dat = STEM_vis.coords_to_grid(input_dat['x'].values,
                                                input_dat['y'].values,
                                                input_dat['COS'].values)
    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(
        os.path.join(input_dir, 'TOPO-124x124.nc'))
    m_pseudo = NAMapFigure(t_str=t_str,
                           cb_axis=cb_axis,
                           map_axis=map_axis)

    m_pseudo.add_ocs_contour_plot(lon,
                                  lat,
                                  gridded_input_dat,
                                  vmin=0.0,
                                  vmax=0.8,
                                  cbar_t_str=cbar_t_str,
                                  colorbar_args={'format': '%0.2f'})
    return(m_pseudo)

def plot_tobspred_conc(input_dir,
                       tobspred_fname,
                       run_dir,
                       map_axis,
                       cb_axis,
                       t_str='',
                       cbar_t_str=''):

    tobspred = STEM_parsers.parse_tobspred(os.path.join(run_dir,
                                                        tobspred_fname))
    gridded_tobspred = STEM_vis.coords_to_grid(tobspred['emi_fac']['x'].values,
                                               tobspred['emi_fac']['y'].values,
                                               tobspred['ocs_conc']['mod'].values)
    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(
        os.path.join(input_dir, 'TOPO-124x124.nc'))
    m_conc = NAMapFigure(t_str=t_str,
                         cb_axis=cb_axis,
                         map_axis=map_axis)

    m_conc.add_ocs_contour_plot(lon,
                                lat,
                                gridded_tobspred,
                                vmin=0.0,
                                vmax=0.8,
                                cbar_t_str=cbar_t_str,
                                colorbar_args={'format': '%0.2f'})
    return(m_conc)

if __name__ == "__main__":
    """
    TODO: make run, input directory arguments
    """

    default_input_dir = '/mnt/home10/skulkarni/StemData21Jul2013/input/'
    default_outfile = os.path.join(os.getenv('HOME'), 'STEM_summary.pdf')

    parser = argparse.ArgumentParser(
        description=("create a one-page PDF with diagnostic plots for " +
                     "a STEM optimization run vs. optimization iteration " +
                     "for a STEM optimization run"))
    parser.add_argument('-r', '--run_dir',
                        nargs='?',
                        type=str,
                        default='',
                        help=('the STEM run directory. Must ' +
                              'contain t_obs_pred*.dat and input.dat'))
    parser.add_argument('-i', '--input_dir',
                        nargs='?',
                        type=str,
                        default=default_input_dir,
                        help=('the STEM input directory. Must ' +
                              'contain TOPO-124x124.nc'))
    parser.add_argument('-o', '--outfile',
                        nargs='?',
                        type=str,
                        default=default_outfile,
                        help='filename for the PDF output file')
    parser.add_argument('-n', '--note',
                        nargs='?',
                        type=str,
                        default='',
                        help=('optional user-specified text to appear' +
                              'in the summary text panel.'))
    args = parser.parse_args()

    matplotlib.rcParams.update({'font.size': 8})
    #plt.rcdefaults # restore defaults


    n_plots = 8
    [fig, ax_list] = initialize_plotting_objects(n_plots)

    print 'plotting summary text panel'
    summary_text_panel(ax_list['text_panel'],
                       args.run_dir,
                       args.input_dir,
                       note=args.note)

    print 'plotting cost function'
    report_opt_df = STEM_parsers.parse_reportopt(os.path.join(args.run_dir,
                                                              'Report.opt'))
    plot_ReportOpt.plot_reportopt(report_opt_df, ax=ax_list['cost_func'])

    print 'plotting emi_fac boxplots'
    emi_fac = STEM_parsers.parse_all_emifac(args.run_dir, mask_ones=False)
    plot_emifac_boxplots.draw_boxplots(emi_fac, ax=ax_list['emi_fac_boxplots'])

    sys.stdout.write('drawing emissions factors map: ')
    file_list = STEM_parsers.get_all_tobspred_fnames(args.run_dir)
    if file_list:
        sys.stdout.write(os.path.basename(file_list[-1]) + '\n')
        sys.stdout.flush()
        file_list = sorted(file_list)
        n_iter = len(file_list)
        midway_through = int(np.floor(n_iter / 2))
        emi_fac_map = plot_tobspred_emifac.draw_plot(
            args.run_dir,
            args.input_dir,
            file_list[-1],
            n_iter,
            t_str='final emissions factors',
            ax=ax_list['emi_fac_map_final'],
            cb_axis=ax_list['emi_fac_map_final_cbar'],
            v_rng=(0.0, 3.0),
            cmap=cm.get_cmap('Oranges'))
        sys.stdout.write('drawing emissions factors map: ' +
                         os.path.basename(file_list[midway_through]) +
                         '\n')
        sys.stdout.flush()
        emi_fac_map = plot_tobspred_emifac.draw_plot(
            args.run_dir,
            args.input_dir,
            file_list[midway_through],
            n_iter,
            t_str='emissions factors, iteration {}'.format(midway_through),
            ax=ax_list['emi_fac_map_N'],
            cb_axis=ax_list['emi_fac_map_N_cbar'],
            v_rng=(0.0, 3.0),
            cmap=cm.get_cmap('Oranges'))

    print 'drawing OCS pseuoddata map'
    psuedodata_OCS_map = plot_inputdat_conc(args.input_dir,
                                            args.run_dir,
                                            ax_list['pseudo_map'],
                                            ax_list['pseudo_map_cbar'],
                                            t_str='pseudodata',
                                            cbar_t_str='OCS [ppbv]')

    print 'drawing OCS prior concentration'
    run_dir_fwd = '/home/thilton/Stem_emi2_onespecies_big_ocssib/run.TWH_fwd_dummy'
    fwd_OCS_map = plot_tobspred_conc(args.input_dir,
                                     't_obs_pred_1.0x.dat',
                                     run_dir_fwd,
                                     ax_list['fwd_map'],
                                     ax_list['fwd_map_cbar'],
                                     t_str='STEM forward run',
                                     cbar_t_str='OCS [ppbv]')

    print 'saving the summary to ' + args.outfile
    fig.savefig(args.outfile)
