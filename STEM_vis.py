"""Provides a number of functions useful for parsing STEM input and
output netcdf files into Python."""

import os.path
import matplotlib.cm as cm
from datetime import datetime, timedelta
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

import na_map
import STEM_parsers

def coords_to_grid(x_seq, y_seq, data_seq, nx=124, ny=124):
    data = np.tile(np.array(np.NaN), [nx, ny])
    for i in range(len(x_seq)):
        data[x_seq[i]-1, y_seq[i]-1] = data_seq[i]

    return(data)

def initialize_STEM_map():
    """Initializes and returns a na_map object with "missing" color
    legend turned off."""
    return(na_map.NAMapFigure(missing_axis=None, cb_axis=True))

def grid_inputdat_data(self, input_dat_file):
    input_dat = STEM_parsers.parse_inputdat(input_dat_file)
    gridded_input_dat = STEM_vis.coords_to_grid(input_dat['x'].values,
                                                input_dat['y'].values,
                                                input_dat['COS'].values)
    return(gridded_input_dat)

def grid_tobspred_data(tobspred):
    """
    Translate parsed columnar t_obs_pred.dat OCS concentrations from a
    data frame to a 2D numpy array.  tobspred is a dict containing two
    pandas data.frames; this will usually be the output of
    STEM_parsers.parse_tobs_pred.

    RETURNS a 2D numpy array of gridded OCS concentrations.
    """
    gridded_tobspred = STEM_vis.coords_to_grid(
        tobspred['emi_fac']['x'].values,
        tobspred['emi_fac']['y'].values,
        tobspred['ocs_conc']['mod'].values)
    return(gridded_tobspred)

def plot_gridded_data(input_dir,
                      gridded_data,
                      map_axis,
                      cb_axis,
                      t_str=' ',
                      cbar_t_str=' ',
                      vmin=0.0,
                      vmax=0.8,
                      extend='neither',
                      cmap=cm.get_cmap('Blues')):
    """
    produce a contour plot on a map for gridded STEM data.  This is
    essentially a wrapper for parsing latitude and longitude and
    collecting the plot labels, colormap, etc.

    RETURNS an na_map object containing the contour plot
    """

    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(
        os.path.join(input_dir, 'TOPO-124x124.nc'))
    map_obj = na_map.NAMapFigure(t_str=t_str,
                                 cb_axis=cb_axis,
                                 map_axis=map_axis)

    map_obj.add_ocs_contour_plot(lon,
                                 lat,
                                 gridded_data,
                                 vmin=vmin,
                                 vmax=vmax,
                                 extend=extend,
                                 cbar_t_str=cbar_t_str,
                                 colorbar_args={'format': '%0.2g'},
                                 cmap=cmap)
    return(map_obj)

def get_midday_mean_ocs_flux(nc_fname):
    """
    calculates mid-day mean OCS surface flux from a I/O API file
    containing OCS flux data.  Mid-day is defined here as 10:00-15:00.

    RETURNS 2D numpy array containing the mean fluxes
    """
    t0 = datetime(2008,6,1)
    ocs_flux = STEM_parsers.parse_STEM_var(nc_fname=nc_fname,
                                           t0=t0,
                                           t1=t0 + timedelta(days=7),
                                           varname='cos')
    idx_midday = np.array(
        [(t.hour >= 10) and (t.hour <= 15) for t in ocs_flux['t']])
    mean_flx = ocs_flux['data'][idx_midday, :, :].mean(axis=0)
    return(mean_flx)
