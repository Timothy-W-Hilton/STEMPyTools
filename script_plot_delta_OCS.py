import os.path
import numpy as np
import pandas as pd
import matplotlib.cm as cm

import STEM_vis
import na_map

ROOT_DIR = os.path.join('/', 'Users', 'Tim', 'work', 'Data',
                        'STEM', 'perturbation_pseudo_data_exp')
# parse input.dat
input_dat_fname = os.path.join( ROOT_DIR, 'input.dat')
input_dat = pd.read_csv(input_dat_fname,
                        sep=None,
                        header=None,
                        skiprows=4,
                        names=('t', 'x', 'y', 'z', 'COS'))

t_obs_pred_1p5X_fname = os.path.join( ROOT_DIR, 't_obs_pred_1.5x.dat')
t_obs_pred_1p5X = pd.read_csv(t_obs_pred_1p5X_fname,
                              delim_whitespace=True,
                              header=None,
                              names=('t', 'obs', 'mod'),
                              nrows=input_dat.shape[0])

t_obs_pred_1p0X_fname = os.path.join( ROOT_DIR, 't_obs_pred_1.0x.dat')
t_obs_pred_1p0X = pd.read_csv(t_obs_pred_1p0X_fname,
                              delim_whitespace=True,
                              header=None,
                              names=('t', 'obs', 'mod'),
                              nrows=input_dat.shape[0])

[lon,lat,topo] =STEM_vis.parse_STEM_coordinates( '/Users/tim/work/Data/STEM/input/TOPO-124x124.nc')

# plot the difference between the two forward model outputs
ocs_1p5 = t_obs_pred_1p5X['mod'].values.reshape(lon.shape)
ocs_1p0 = t_obs_pred_1p0X['mod'].values.reshape(lon.shape)

m_1p5 = na_map.NAMapFigure(t_str="[OCS], 1.5x CASA sfc flux",
                       missing_axis=False)
cs = m_1p5.add_ocs_contour_plot(lon,
                               lat,
                               ocs_1p0,
                               cbar_t_str='[OCS], ppbv')

m_1p0 = na_map.NAMapFigure(t_str="[OCS], 1.0x CASA sfc flux",
                       missing_axis=False)
cs_1p0 = m_1p0.add_ocs_contour_plot(lons=lon,
                               lats=lat,
                               data=ocs_1p0,
                               cbar_t_str='[OCS], ppbv')

delta_ocs = ocs_1p5 - ocs_1p0
m_diff = na_map.NAMapFigure(t_str="$\Delta$[OCS], 1.5x - 1.0x CASA sfc flux",
                           missing_axis=False)
cmap_diff = cm.get_cmap('winter')
cs_diff = m_diff.add_ocs_contour_plot(lons=lon,
                                     lats=lat,
                                     data=delta_ocs,
                                     cmap=cmap_diff,
                                     vmax=np.abs(delta_ocs).max(),
                                     vmin=-np.abs(delta_ocs).max(),
                                     cbar_t_str='[OCS], ppbv')

m_gt0 = na_map.NAMapFigure(t_str="$\Delta$[OCS] > 0",
                           missing_axis=False)
delta_gt0 = np.ma.masked_array(delta_ocs, delta_ocs < 0)
cs_diff = m_gt0.add_ocs_contour_plot(lons=lon,
                                     lats=lat,
                                     data=delta_gt0.mask,
                                     cbar_t_str='$\Delta$[OCS] > 0')

# make sure the STEM grid is oriented correctly -- the TOPO input
# netcdf file says in its attributes that the origin is in the SE
# corner.  This X should increase E to W and y should increase S to N.
x = input_dat['x'].values.reshape(lon.shape)
y = input_dat['y'].values.reshape(lon.shape)
m_x = na_map.NAMapFigure(t_str="STEM grid x values",
                           missing_axis=False)
cs_x = m_x.add_ocs_contour_plot(lons=lon,
                                lats=lat,
                                data=x,
                                cmap=cm.get_cmap('binary'),
                                cbar_t_str="X")
m_y = na_map.NAMapFigure(t_str="STEM grid y values",
                           missing_axis=False)
cs_y = m_y.add_ocs_contour_plot(lons=lon,
                                lats=lat,
                                data=y,
                                cmap=cm.get_cmap('binary'),
                                cbar_t_str="Y")
