"""Provides a number of functions useful for parsing STEM input and
output netcdf files into Python."""

from datetime import datetime
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

import na_map

def coords_to_grid(x_seq, y_seq, data_seq, nx=124, ny=124):
    data = np.tile(np.array(np.NaN), [nx, ny])
    for i in range(len(x_seq)):
        data[x_seq[i]-1, y_seq[i]-1] = data_seq[i]

    return(data)

def initialize_STEM_map():
    """Initializes and returns a na_map object with "missing" color
    legend turned off."""
    return(na_map.NAMapFigure(missing_axis=None, cb_axis=True))
