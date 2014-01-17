"""
short module to produce contour plots of STEM grid X values and STEM
grid Y values.  These are useful to make sure that data are being read
and plotted with the correct orientation.  The minimum X and Y values
should be in the southeast corner of the domain.
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt

import STEM_vis
import na_map

if __name__ == "__main__":

    INPUT_DIR = os.path.join('/mnt',
                             'home10',
                             'skulkarni',
                             'StemData21Jul2013',
                             'input')
    nx = 124
    ny = 124

    xvals = np.repeat(np.arange(nx), ny)
    yvals = np.tile(np.arange(ny), nx)

    xgrid = STEM_vis.coords_to_grid(xvals, yvals, xvals)
    ygrid = STEM_vis.coords_to_grid(xvals, yvals, yvals)

    [lon,lat,topo] =STEM_vis.parse_STEM_coordinates(
        os.path.join(INPUT_DIR, 'TOPO-124x124.nc'))

    m_x = na_map.NAMapFigure(t_str="STEM 124x124 X grid",
                             cb_axis=True)
    m_x.add_ocs_contour_plot(lon, lat, xgrid, colorbar_args={'format':'%d'})

    m_y = na_map.NAMapFigure(t_str="STEM 124x124 Y grid",
                             cb_axis=True)
    m_y.add_ocs_contour_plot(lon, lat, ygrid, colorbar_args={'format':'%d'})
