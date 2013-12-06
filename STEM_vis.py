"""Provides a number of functions useful for parsing STEM input and
output netcdf files into Python."""

from datetime import datetime
import numpy as np
from netCDF4 import Dataset

import na_map

def parse_STEM_var(nc_fname=None, t0=None, t1=None, varname=None):
    """ Parse a STEM variable from a STEM I/O API netcdf file.
    varname must be a variable in the netcdf file. The file must also
    contain a variable TFLAG containing timestamps in the format
    <YYYYDDD,HHMMSS>.  Returns the values in varname as well as the
    timestamps (as datenum.datenum objects)."""
    nc = Dataset(nc_fname, 'r', format='NETCDF4')
    # read timestamps to datetime.datetime
    t = np.squeeze(nc.variables['TFLAG'])
    t_dt = np.array(([datetime.strptime(str(this[0]) +
                                        str(this[1]).zfill(6), '%Y%j%H%M%S')
                                        for this in t]))
    # find the requested timestamps
    t_idx = (t_dt >= t0) & (t_dt <= t1)
    # retrieve the requested [OCS] data
    data = nc.variables[varname][t_idx, 0, :, : ]
    return( data, t_dt[t_idx] )

def parse_STEM_coordinates(topo_fname):
    """Parse STEM grid latitude and longitude."""
    topo = Dataset(topo_fname, 'r', format='NETCDF4')

    lat = np.squeeze(topo.variables['LAT'])
    lon = np.squeeze(topo.variables['LON'])

    return(lon, lat)

def initialize_STEM_map():
    """Initializes and returns a na_map object with "missing" color
    legend turned off."""
    return(na_map.NAMapFigure(missing_axis=False))


