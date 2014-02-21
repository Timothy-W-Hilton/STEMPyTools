"""
produce a one-page PDF showing some plots that describe the twin
experiments setup
"""

import matplotlib
matplotlib.use('Agg') # this seems to get around the 'couldn't connect
                        # to display "localhost:10.0" error' that
                        # crops up now and then -- TWH
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
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

class twinExperimentSummary(object):
    """class to provide a container for stuff needed for summary
    figure of twin experiment"""

    def __init__(self,
                 outfile_fname,
                 input_dir):

        self.input_dir = input_dir
        self.outfile_fname = outfile_fname

        self.fig, self.axes_list = self.initialize_plotting_objects()

        #specify the files containing the OCS concentrations for the
        #1.0x, 1.5x, and 0.8x run
        self.ocs_conc_1p0x = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'run.TWH_fwd_dummy',
            't_obs_pred_1.0x.dat')
        self.ocs_conc_1p5x = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'run.TWH_opt_test_large_slab_strong_prior',
            'input.dat')
        self.ocs_conc_0p8x = os.path.join(
            os.getenv('HOME'),
            'Stem_emi2_onespecies_big_ocssib',
            'run.TWH_opt_test_large_slab_weak_prior_0.8',
            'input.dat')

        #get the gridded concentrations
        self.gridded_ocs_1p0x = self.grid_tobspred_data(self.ocs_conc_1p0x)
        self.gridded_ocs_1p5x = self.grid_inputdat_data(self.ocs_conc_1p5x)
        self.gridded_ocs_0p8x = self.grid_inputdat_data(self.ocs_conc_0p8x)

        all_ocs = np.concatenate([self.gridded_ocs_1p0x,
                                 self.gridded_ocs_1p5x,
                                 self.gridded_ocs_0p8x])
        ocs_rng = (all_ocs.min(), all_ocs.max())

        # plot the gridded concentrations
        print 'drawing [OCS] forward run map'
        self.plot_gridded_data(self.input_dir,
                               self.gridded_ocs_1p0x,
                               self.axes_list['ocs_conc_1.0_map'],
                               self.axes_list['ocs_conc_1.0_cbar'],
                               t_str='A. STEM forward run\n("true" [OCS])',
                               cbar_t_str='[OCS] (ppbv)',
                               vmin=ocs_rng[0],
                               vmax=ocs_rng[1])
        print 'drawing [OCS] forward run 1.5x map'
        self.plot_gridded_data(self.input_dir,
                               self.gridded_ocs_1p5x,
                               self.axes_list['ocs_conc_1.5_map'],
                               self.axes_list['ocs_conc_1.5_cbar'],
                               t_str=('B. pseudo-obervations, F = 1.5\n' +
                                      '(prior [OCS], experiments 1,2,3)'),
                               cbar_t_str='[OCS] (ppbv)',
                               vmin=ocs_rng[0],
                               vmax=ocs_rng[1])
        print 'drawing [OCS] forward run 0.8x map'
        self.plot_gridded_data(self.input_dir,
                               self.gridded_ocs_0p8x,
                               self.axes_list['ocs_conc_0.8_map'],
                               self.axes_list['ocs_conc_0.8_cbar'],
                               t_str=('D. pseudo-observations, F = 0.8\n' +
                                      '(prior [OCS], experiment 4)'),
                               cbar_t_str='[OCS] (ppbv)',
                               vmin=ocs_rng[0],
                               vmax=ocs_rng[1])
        print 'drawing OCS flux prior map'
        prior_flx = STEM_vis.get_midday_mean_ocs_flux(
            nc_fname=os.path.join(self.input_dir,
                                  'surfem-124x124-casa-cos_2008_2009.nc'))
        self.plot_gridded_data(
            self.input_dir,
            prior_flx,
            self.axes_list['ocs_prior_flux_map'],
            self.axes_list['ocs_prior_flux_cbar'],
            t_str=('F. prior OCS flux (from CASA)'),
            cbar_t_str=r'mol OCS m$^{-2}$ s$^{-1}$',
            vmin=-2.5e-10,
            vmax=prior_flx.max(),
            cmap=cm.get_cmap('Greens_r'),
            colorbar_args={'format':'%0.1e'})
        # print 'drawing OCS flux posterior map'
        # self.plot_gridded_data(
        #     self.input_dir,
        #     self.get_midday_mean_ocs_flux(
        #         nc_fname=os.path.join(
        #             self.run_dir,
        #             output,
        #             'AQOUT-124x124-22levs-casa-cos_2008_2009.nc')),
        #     self.axes_list['ocs_post_flux_map'],
        #     self.axes_list['ocs_post_flux_cbar'],
        #     t_str=('G. OCS flux posterior'),
        #     cbar_t_str='OCS flux')

        #plot delta OCS contours
        delta_1p5_1p0 = self.gridded_ocs_1p5x - self.gridded_ocs_1p0x
        delta_1p0_0p8 = self.gridded_ocs_1p0x - self.gridded_ocs_0p8x
        all_delta = np.concatenate((delta_1p0_0p8, delta_1p5_1p0))
        delta_rng = (all_delta.min(), all_delta.max())
        print 'drawing delta [OCS], 1.5x - 1.0x'
        self.plot_gridded_data(
            self.input_dir,
            delta_1p5_1p0,
            self.axes_list['ocs_conc_1.5_1.0_delta_map'],
            self.axes_list['ocs_conc_1.5_1.0_delta_cbar'],
            t_str=r'C. $\Delta$ [OCS], F=1.5 - F=1.0 ',
            cbar_t_str=r'$\Delta$ [OCS] (ppbv)',
            vmin=delta_rng[0],
            vmax=delta_rng[1],
            cmap=cm.get_cmap('Reds'))
        print 'drawing delta [OCS], 1.0x - 0.8x'
        self.plot_gridded_data(
            self.input_dir,
            delta_1p0_0p8,
            self.axes_list['ocs_conc_1.0_0.8_delta_map'],
            self.axes_list['ocs_conc_1.0_0.8_delta_cbar'],
            t_str=r'E. $\Delta$ [OCS], F=1.0 - F=0.8x',
            cbar_t_str=r'$\Delta$ [OCS] (ppbv)',
            vmin=delta_rng[0],
            vmax=delta_rng[1],
            cmap=cm.get_cmap('Reds'))

        print 'saving the summary to ' + args.outfile
        self.fig.savefig(self.outfile_fname)

    def initialize_plotting_objects(self, n_plots=8, figsize=(8.5, 11)):
        """create a figure and axes instances to hold STEM run summary
        plots.  The axes are arranged in two columns, with the number of
        rows determined by the number of plots requested. Returns a tuple
        containing the figure and a list of axes"""

        n_rows = np.int(np.ceil(n_plots/2))
        #This is a little bit of gridspec finagling -- there will be two
        #columns of plots in the final product, but specifying 10 allows
        #each column to be subdivided into five "regions" (two columns
        #times five regions equals ten).  Then within each subdivided
        #column the the first four regions will contain the map, and the
        #fifth the colorbar.
        n_cols = 10

        fig = plt.figure(figsize=figsize)
        #define spacing for left-hand column of panels

        gs1 = gridspec.GridSpec(n_rows, n_cols)
        gs1.update(left=0.10, right=0.50, wspace=0.05, hspace=0.30)
        ax_list = {'ocs_conc_1.0_map':fig.add_subplot(gs1[0, 0:7]),
                   'ocs_conc_1.0_cbar':fig.add_subplot(gs1[0, 8]),
                   'ocs_conc_1.5_map':fig.add_subplot(gs1[1, 0:7]),
                   'ocs_conc_1.5_cbar':fig.add_subplot(gs1[1, 8]),
                   'ocs_conc_0.8_map':fig.add_subplot(gs1[2, 0:7]),
                   'ocs_conc_0.8_cbar':fig.add_subplot(gs1[2, 8]),
                   'ocs_prior_flux_map':fig.add_subplot(gs1[3, 0:7]),
                   'ocs_prior_flux_cbar':fig.add_subplot(gs1[3, 8])}
        #define spacing for right-hand column of panels
        gs2 = gridspec.GridSpec(n_rows, n_cols)
        gs2.update(left=0.55, right=0.95, wspace=0.05, hspace=0.30)
        ax_list.update({'ocs_conc_1.5_1.0_delta_map':fig.add_subplot(gs2[1, 0:7]),
                        'ocs_conc_1.5_1.0_delta_cbar':fig.add_subplot(gs2[1, 8]),
                        'ocs_conc_1.0_0.8_delta_map':fig.add_subplot(gs2[2, 0:7]),
                        'ocs_conc_1.0_0.8_delta_cbar':fig.add_subplot(gs2[2, 8])})
                        # 'ocs_post_flux_map':fig.add_subplot(gs2[3, 0:7]),
                        # 'ocs_post_flux_cbar':fig.add_subplot(gs2[3, 8])})

        return((fig, ax_list))

    def grid_inputdat_data(self, input_dat_file):
        input_dat = STEM_parsers.parse_inputdat(input_dat_file)
        gridded_input_dat = STEM_vis.coords_to_grid(input_dat['x'].values,
                                                    input_dat['y'].values,
                                                    input_dat['COS'].values)
        return(gridded_input_dat)

    def grid_tobspred_data(self, tobspred_fname):
        tobspred = STEM_parsers.parse_tobspred(tobspred_fname)
        gridded_tobspred = STEM_vis.coords_to_grid(
            tobspred['emi_fac']['x'].values,
            tobspred['emi_fac']['y'].values,
            tobspred['ocs_conc']['mod'].values)
        return(gridded_tobspred)

    def plot_gridded_data(self,
                          input_dir,
                          gridded_data,
                          map_axis,
                          cb_axis,
                          t_str='',
                          cbar_t_str='',
                          vmin=0.0,
                          vmax=0.8,
                          cmap=cm.get_cmap('Blues'),
                          colorbar_args={'format': '%0.2f'}):

        lon, lat, topo = STEM_parsers.parse_STEM_coordinates(
            os.path.join(input_dir, 'TOPO-124x124.nc'))
        map_obj = NAMapFigure(t_str=t_str,
                              cb_axis=cb_axis,
                              map_axis=map_axis)

        map_obj.add_ocs_contour_plot(lon,
                                     lat,
                                     gridded_data,
                                     vmin=vmin,
                                     vmax=vmax,
                                     cbar_t_str=cbar_t_str,
                                     colorbar_args=colorbar_args,
                                     cmap=cmap)
        return(map_obj)

if __name__ == "__main__":

    default_input_dir = '/mnt/home10/skulkarni/StemData21Jul2013/input/'
    default_outfile = os.path.join(os.getenv('PLOTS'),
                                   'STEM_twin_exp_summary.pdf')
    #handle input arguments, provide online documentation
    parser = argparse.ArgumentParser(
        description=("create a one-page PDF with some plots that " +
                     "summarize the twin experiment setup"))
    parser.add_argument('-i', '--input_dir',
                        nargs='?',
                        type=str,
                        default=default_input_dir,
                        help=('the STEM input directory. Must ' +
                              'contain TOPO-124x124.nc. Default is ' +
                              default_input_dir))
    parser.add_argument('-o', '--outfile',
                        nargs='?',
                        type=str,
                        default=default_outfile,
                        help=('filename for the PDF output file. ' +
                              'Default is ' + default_outfile))
    args = parser.parse_args()

    matplotlib.rcParams.update({'font.size': 8})
    summary_obj = twinExperimentSummary(outfile_fname=args.outfile,
                                        input_dir=args.input_dir)
