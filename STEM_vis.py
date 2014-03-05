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

def grid_inputdat_data(input_dat_data):
    """
    place input.dat OCS concentrations into a 2D numpy array.
    input_dat_data is the output of STEM_parsers.parse_inputdat().
    """
    gridded_input_dat = coords_to_grid(input_dat_data['x'].values,
                                       input_dat_data['y'].values,
                                       input_dat_data['COS'].values)
    return(gridded_input_dat)

def grid_tobspred_data(tobspred, which_data='ocs_mod'):
    """
    Translate parsed columnar t_obs_pred.dat OCS concentrations from a
    data frame to a 2D numpy array.  tobspred is a dict containing two
    pandas data.frames; this will usually be the output of
    STEM_parsers.parse_tobs_pred.

    which_data: {'ocs_mod'}|'ocs_obs'|'emi_fac'; allows the caller to
    request either OCS model concentrations, OCS observed
    concentrations, or emissions factors from the t_obs_pred.dat file.

    RETURNS a 2D numpy array of gridded OCS concentrations.
    """

    if which_data is 'ocs_obs':
        data = tobspred['ocs_conc']['obs'].values
    elif which_data is 'ocs_mod':
        data = tobspred['ocs_conc']['mod'].values
    elif which_data is 'emi_fac':
        data = tobspred['emi_fac']['emi_fac'].values

    gridded_tobspred = coords_to_grid(tobspred['emi_fac']['x'].values,
                                      tobspred['emi_fac']['y'].values,
                                      data)
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

def calc_gridded_STEM_errors(top_data):
    """
    calculate gridded background and mismatch error from a STEM
    t_obs_pred*.dat file.  Background error is defined as (calculated
    emission factor minus 1.0).  Mismatch error is defined as (model
    [OCS] minus observed [OCS]).

    'top_data' abbreviates 't_obs_pred data'.

    RETURNS a dict with keys 'mismatch' and 'background'; each refer
    to a 2D numpy array of gridded errors.
    """

    gr_mod = grid_tobspred_data(top_data, 'ocs_mod')
    gr_obs = grid_tobspred_data(top_data, 'ocs_obs')
    gr_mismatch = gr_mod - gr_obs

    gr_emi = grid_tobspred_data(top_data, 'emi_fac')
    gr_background = gr_emi - np.ones_like(gr_emi)

    return({'mismatch':gr_emi, 'background':gr_background})
