### produce a contour plot of scaling factors from a STEM t_obs_pred.dat file

import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from STEM_parsers import parse_inputdat, parse_tobspred
import STEM_vis
import na_map

if __name__ == "__main__":

    RUN_DIR = os.path.join('/home',
                            'thilton',
                            'Stem_emi2_onespecies_big_ocssib',
                            'run.TWH_opt_test_large_slab')
    INPUT_DIR = os.path.join('/mnt',
                             'home10',
                             'skulkarni',
                             'StemData21Jul2013',
                             'input')

    # RUN_DIR = os.path.join('/', 'Users', 'Tim', 'work', 'Data',
    #                         'STEM', 'perturbation_pseudo_data_exp')
    #parse input.dat
    inputdat_fname = os.path.join( RUN_DIR, 'input.dat')
    inputdat = parse_inputdat(inputdat_fname)
    #parse t_obs_pred.dat emissions factors
    tobspred_fname = os.path.join( RUN_DIR, 't_obs_pred.dat')
    tobspred = parse_tobspred(tobspred_fname)['emi_fac']

    # translate t_obs_pred emi_fac values into 124 x 124 grid
    [lon,lat,topo] =STEM_vis.parse_STEM_coordinates(
        os.path.join(INPUT_DIR, 'TOPO-124x124.nc'))
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
