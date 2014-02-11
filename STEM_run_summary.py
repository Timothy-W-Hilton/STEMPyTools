import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import os.path
import glob

import STEM_parsers
import plot_ReportOpt
import plot_emifac_boxplots
import plot_tobspred_emifac

def initialize_plotting_objects(n_plots=8, figsize=(17,22)):
    """create a figure and axes instances to hold STEM run summary
    plots.  Returns a tuple containing the figure and a list of axes"""

    n_rows = np.int(np.ceil(n_plots/2))

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(n_rows, 2)
    #ax = [plt.subplot(gs[i]) for i in range(n_plots)]
    ax_list = [fig.add_subplot(gs[i]) for i in range(n_plots)]

    return((fig, ax_list))


if __name__ == "__main__":

    run_dir = os.path.join(os.getenv('HOME'),
                           'Stem_emi2_onespecies_big_ocssib',
                           'run.TWH_opt_test_large_slab_weak_prior')
    input_dir = '/mnt/home10/skulkarni/StemData21Jul2013/input/'

    n_plots = 8
    [fig, ax_list] = initialize_plotting_objects(n_plots)

    print 'plotting cost function'
    report_opt_df = STEM_parsers.parse_reportopt(os.path.join(run_dir,
                                                              'Report.opt'))
    plot_ReportOpt.plot_reportopt(report_opt_df, ax=ax_list[0])

    print 'plotting emi_fac boxplots'
    emi_fac = STEM_parsers.parse_all_emifac(run_dir, mask_ones=False)
    plot_emifac_boxplots.draw_boxplots(emi_fac, ax=ax_list[1])

    print 'drawing emissions factors map'
    file_list = glob.glob(os.path.join(run_dir, 't_obs_pred*.dat'))
    if file_list:
        file_list = sorted(file_list)
        n_iter = len(file_list)
        emi_fac_map = plot_tobspred_emifac.draw_plot(run_dir,
                                                     input_dir,
                                                     file_list[-1],
                                                     n_iter,
                                                     ax=ax_list[2],
                                                     cb_axis=ax_list[3],
                                                     v_rng=(0.0, 2.0))
    print 'saving the summary'
    fig.savefig(os.path.join(os.getenv('HOME'), 'Plots', 'test_summary.pdf'))
    
