import numpy as np
import os
import os.path

import STEM_parsers
import STEM_vis

def get_conc_ratios_vs_fwd(run_dir):

    input_dat = STEM_parsers.parse_inputdat(
        os.path.join(os.getenv('HOME'),
                     'Stem_emi2_onespecies_big_ocssib',
                     'run.TWH_opt_test_large_slab_weak_prior_0.8',
                     'input.dat'))
    gridded_ptbd = STEM_vis.coords_to_grid(input_dat['x'].values,
                                           input_dat['y'].values,
                                           input_dat['COS'].values)

    # ptbd_conc = STEM_parsers.parse_tobspred(
    #     os.path.join(os.getenv('HOME'),
    #                  'Stem_emi2_onespecies_big_ocssib',
    #                  'run.TWH_opt_test_large_slab_13Feb',
    #                  'input.dat'))
    # gridded_ptbd = STEM_vis.coords_to_grid(ptbd_conc['emi_fac']['x'].values,
    #                                            ptbd_conc['emi_fac']['y'].values,
    #                                            ptbd_conc['ocs_conc']['mod'].values)

    fwd_conc = STEM_parsers.parse_tobspred(
        os.path.join(os.getenv('HOME'),
                     'Stem_emi2_onespecies_big_ocssib',
                     'run.TWH_fwd_dummy',
                     't_obs_pred_1.0x.dat'))
    gridded_fwd = STEM_vis.coords_to_grid(fwd_conc['emi_fac']['x'].values,
                                                   fwd_conc['emi_fac']['y'].values,
                                                   fwd_conc['ocs_conc']['mod'].values)

    ratio = gridded_ptbd / gridded_fwd
    delta = gridded_ptbd - gridded_fwd

    return({'ratio':ratio, 'delta':delta})
