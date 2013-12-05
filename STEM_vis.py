"""Provides a number of functions useful for parsing STEM input and
output netcdf files into Python."""

from datetime import datetime
import numpy as np
from netCDF4 import Dataset

import na_map

def parse_STEM_AQOUT(aq_fname=None, t0=None, t1=None):
    """ Parse a OCS from a STEM AQOUT file.  Assumes that the OCS
    variable in the netcdf file is called CO2_TRACER1.  Returns OCS
    and timestamps (as datenum.datenum objects)."""
    aq_out = Dataset(aq_fname, 'r', format='NETCDF4')
    # read timestamps to datetime.datetime
    t = np.squeeze(aq_out.variables['TFLAG'])
    t_dt = np.array(([datetime.strptime(str(this[0]) +
                                        str(this[1]).zfill(6), '%Y%j%H%M%S')
                                        for this in t]))
    # find the requested timestamps
    t_idx = (t_dt >= t0) & (t_dt <= t1)
    # retrieve the requested [OCS] data
    ocs = aq_out.variables['CO2_TRACER1'][t_idx, 0, :, : ]
    return( ocs, t_dt[t_idx] )

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


