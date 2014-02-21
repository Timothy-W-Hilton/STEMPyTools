import matplotlib
matplotlib.use('Agg')   # this seems to get around the 'couldn't connect to display "localhost:10.0" error' that crops up now and then -- TWH

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib import rc
import numpy as np
import os
import os.path
import pdb
import sys
from datetime import datetime, timedelta
import argparse

import STEM_parsers
import STEM_vis
from na_map import NAMapFigure
import plot_ReportOpt
import plot_emifac_boxplots
import plot_tobspred_emifac

class STEMRun(object):
    """
    container class for data useful for summarizing a STEM run.
    Attributes are:
    - run_dir: str; full path to the run directory for the STEM run
    - input_dir: str; full path to the directory containing the run's
      input files
    - prior_flux_fname: str; full path to the input file specifying
      prior OCS flux
    - coord_fname: str; full path to the input file specifying STEM
      grid latitude & longitude
    - top_fnames: str; full paths of all t_obs_pred*.dat files from the STEM run
    - prior_ocs_flux: 2D numpy array; prior OCS flux [mol m-2 s-1]
    - post_ocs_flux: 2D numpy array; posterior OCS flux [mol m-2 s-1]
    - prior_ocs_conc: 2D numpy array; prior OCS concentration [ppbv]
    - post_ocs_conc: 2D numpy array; posterior OCS concentration [ppbv]

    """

    def __init__(self,
                 input_dir,
                 run_dir,
                 prior_flx_fname='surfem-124x124-casa-cos_2008_2009.nc',
                 coord_fname = 'TOPO-124x124.nc'):
        '''
        populate a STEMRun object from a specified run directory and
        input directory.  prior_flux_fname specifies the basename of
        the I/O API file containing prior OCS fluxes.  coord_fname
        specifies the basename of the I/O API file containing the
        latitudes and longitudes of the STEM grid cell centers.
        '''
        self.input_dir = input_dir
        self.run_dir = run_dir
        self.prior_flx_fname=os.path.join(self.input_dir,
                                          prior_flx_fname)
        self.coord_fname=os.path.join(self.input_dir,
                                      coord_fname)

        # get intial and final t_obs_pred.dat file names
        self.top_fnames = STEM_parsers.get_all_tobspred_fnames(
            self.run_dir)

        fnl_top_data = STEM_parsers.parse_tobspred(self.top_fnames[-1])
        self.final_emifac = STEM_vis.grid_tobspred_data(fnl_top_data,
                                                        which_data='emi_fac')

        self.prior_ocs_flux = STEM_vis.get_midday_mean_ocs_flux(
            self.prior_flx_fname)
        self.post_ocs_flux = self.prior_ocs_flux * self.final_emifac

        ini_top_data = STEM_parsers.parse_tobspred(self.top_fnames[0])
        self.prior_ocs_conc = STEM_vis.grid_tobspred_data(ini_top_data,
                                                          which_data='ocs_obs')
        self.post_ocs_conc = STEM_vis.grid_tobspred_data(fnl_top_data,
                                                         which_data='ocs_mod')


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
                    'post_conc_map':fig.add_subplot(gs1[3, 0:7]),
                    'post_conc_cbar':fig.add_subplot(gs1[3, 8])})
    #define spacing for right-hand column of panels
    gs2 = gridspec.GridSpec(n_rows, n_cols)
    gs2.update(left=0.55, right=0.95, wspace=0.05, hspace=0.30)
    ax_list.update({'emi_fac_boxplots':fig.add_subplot(gs2[1, 0:9]),
                    'emi_fac_map_N':fig.add_subplot(gs2[2, 0:7]),
                    'emi_fac_map_N_cbar':fig.add_subplot(gs2[2, 8]),
                    'post_flux_map':fig.add_subplot(gs2[3, 0:7]),
                    'post_flux_cbar':fig.add_subplot(gs2[3, 8])})

    return((fig, ax_list))

def summary_text_panel(ax, run_dir, input_dir, fontsize=8, note=' '):
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
                       t_str=' ',
                       cbar_t_str=' '):

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
                       t_str=' ',
                       cbar_t_str=' '):

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
                                vmin=gridded_tobspred.min(),
                                vmax=gridded_tobspred.max(),
                                cbar_t_str=cbar_t_str,
                                colorbar_args={'format': '%0.2f'})
    return(m_conc)

def get_AQOUT_midday_mean_ocs_flux(AQOUT_fname):
    """
    Parses an OCS flux I/O API file (e.g. CASA surface fluxes or a
    STEM AQOUT file) and caculates the mean OCS fluxes for the hours
    10:00 to 15:00.  The input file must contain the netcdf variables T_FLAG and ocs

    RETURNS a 2D numpy array of gridded mean OCS fluxes.
    """
    t0 = datetime(2008,6,1)
    ocs_flux = STEM_parsers.parse_STEM_var(nc_fname=AQOUT_fname,
                                           t0=t0,
                                           t1=t0 + timedelta(days=7),
                                           varname='CO2_TRACER1')
    idx_midday = np.array(
        [(t.hour >= 10) and (t.hour <= 15) for t in ocs_flux['t']])
    mean_flx = ocs_flux['data'][idx_midday, :, :].mean(axis=0)
    return(mean_flx)

def create_summary(this_run,
                   note,
                   outfile='STEM_run_summary.pdf',
                   n_plots=8):
    [fig, ax_list] = initialize_plotting_objects(n_plots)

    print 'plotting summary text panel'
    summary_text_panel(ax_list['text_panel'],
                       this_run.run_dir,
                       this_run.input_dir,
                       note=note)

    print 'plotting cost function'
    report_opt_df = STEM_parsers.parse_reportopt(os.path.join(this_run.run_dir,
                                                              'Report.opt'))
    plot_ReportOpt.plot_reportopt(report_opt_df, ax=ax_list['cost_func'])

    print 'plotting emi_fac boxplots'
    emi_fac = STEM_parsers.parse_all_emifac(this_run.run_dir, mask_ones=False)
    plot_emifac_boxplots.draw_boxplots(emi_fac, ax=ax_list['emi_fac_boxplots'])

    sys.stdout.write('drawing posterior emissions factors map:')
    if this_run.top_fnames:
        sys.stdout.write(os.path.basename(this_run.top_fnames[-1]) + '\n')
        sys.stdout.flush()
        n_iter = len(this_run.top_fnames)
        emi_fac_map = plot_tobspred_emifac.draw_plot(
            this_run.run_dir,
            this_run.input_dir,
            this_run.top_fnames[-1],
            n_iter,
            t_str='final emissions factors',
            ax=ax_list['emi_fac_map_final'],
            cb_axis=ax_list['emi_fac_map_final_cbar'],
            v_rng=(0.0, None),
            extend='max',
            cmap=cm.get_cmap('Oranges'))


    print 'drawing OCS posterior concentration'
    run_dir_fwd = '/home/thilton/Stem_emi2_onespecies_big_ocssib/run.TWH_fwd_dummy'
    final_tobs_pred = STEM_parsers.get_all_tobspred_fnames(this_run.run_dir)[-1]
    post_ocs_conc_map = STEM_vis.plot_gridded_data(this_run.input_dir,
                                                   this_run.post_ocs_conc,
                                                   ax_list['post_conc_map'],
                                                   ax_list['post_conc_cbar'],
                                                   t_str='posterior [OCS]',
                                                   cbar_t_str='OCS [ppbv]')
    print 'drawing OCS posterior flux'
    post_ocs_flux_map = STEM_vis.plot_gridded_data(
        this_run.input_dir,
        this_run.post_ocs_flux,
        ax_list['post_flux_map'],
        ax_list['post_flux_cbar'],
        t_str='posterior OCS flux',
        cbar_t_str=r'mol OCS m$^{-2}$ s$^{-1}$',
        cmap=cm.get_cmap('Greens_r'),
        vmin=-2.5e-10,
        vmax=0.0,
        extend='min')

    print 'saving the summary to ' + outfile
    fig.savefig(args.outfile)


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
                        default=' ',
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
                        default=' ',
                        help=('optional user-specified text to appear' +
                              'in the summary text panel.'))
    args = parser.parse_args()

    this_run = STEMRun(args.input_dir,
                       args.run_dir)

    matplotlib.rcParams.update({'font.size': 8})
    #plt.rcdefaults # restore defaults

    create_summary(this_run,
                   args.note,
                   args.outfile)
