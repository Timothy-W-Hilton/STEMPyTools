'''adapted from http://matplotlib.org/examples/pylab_examples/multipage_pdf.html'''

import os.path
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse
import glob
import re

from plot_tobspred_emifac import draw_plot

if __name__ == '__main__':

    # define default directories, files
    default_input_dir = os.path.join('/mnt',
                                     'home10',
                                     'skulkarni',
                                     'StemData21Jul2013',
                                     'input')
    default_run_dir = os.path.join('/home',
                                   'thilton',
                                   'Stem_emi2_onespecies_big_ocssib',
                                   'run.TWH_opt_test_large_slab_weak_prior')
    default_fname_pdf = os.path.join('/home',
                                     'thilton',
                                     'Plots',
                                     'large_slab_weakpriors_emifac_map.pdf')

    # setup an argument parser
    parser = argparse.ArgumentParser(
        description=("draw a contour plots of scaling factors from " +
                     "STEM t_obs_pred.dat files from all iterations of " +
                     "an optimization run.  The plots are sent to a " +
                     "multi-page PDF file."))
    parser.add_argument('-r', '--run_dir',
                        nargs='?',
                        type=str,
                        default=os.getcwd(),
                        help=('the STEM run directory. Must ' +
                              'contain at least one t_obs_pred_ABC.dat, ' +
                              'with ABC a three-digit integer indicating' +
                              'iteration.'))
    parser.add_argument('-i', '--input_dir',
                        nargs='?',
                        type=str,
                        default=default_input_dir,
                        help=('the STEM input directory. Must ' +
                              'contain TOPO-124x124.nc'))
    parser.add_argument('-o', '--outfilename',
                        nargs='?',
                        type=str,
                        default='./STEM_emifac_contour_maps.pdf',
                        dest='outfilename',
                        help=('the name of the pdf file containing ' +
                              'all output plots'))
    args = parser.parse_args()

    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    pdf = PdfPages(args.outfilename)
    file_list = glob.glob(os.path.join(args.run_dir,
                                       't_obs_pred*.dat'))
    for this_file in file_list:
        match_obj = re.search('[0-9]{3}', this_file)
        if match_obj:
            n_iter = int(match_obj.group())
            print 'drawing countours for ' + this_file
            m = draw_plot(args.run_dir,
                          args.input_dir,
                          this_file,
                          n_iter)
            pdf.savefig(m.fig)
            plt.close(m.fig)

    pdf.close()
