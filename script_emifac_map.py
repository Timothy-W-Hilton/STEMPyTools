'''adapted from http://matplotlib.org/examples/pylab_examples/multipage_pdf.html'''

import os.path
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from plot_tobspred_emifac import draw_plot

input_dir = os.path.join('/mnt',
                         'home10',
                         'skulkarni',
                         'StemData21Jul2013',
                         'input')
run_dir = os.path.join('/home',
                         'thilton',
                         'Stem_emi2_onespecies_big_ocssib',
                         'run.TWH_opt_test_large_slab_weak_prior')

fname_pdf = os.path.join('/home',
                         'thilton',
                         'Plots',
                         ('large_slab_weakpriors_emifac_map.pdf'))

# Create the PdfPages object to which we will save the pages:
# The with statement makes sure that the PdfPages object is closed properly at
# the end of the block, even if an Exception occurs.
pdf = PdfPages(fname_pdf)
for this_iter in range(1,20):
    m = draw_plot(run_dir,
                  input_dir,
                  't_obs_pred_{:03d}.dat'.format(this_iter),
                  this_iter)
    print 'saving iteration {}'.format(this_iter) 
    pdf.savefig(m.fig)
    plt.close(m.fig)

pdf.close()
