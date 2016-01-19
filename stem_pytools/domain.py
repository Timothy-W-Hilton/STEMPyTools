import os
import os.path
import matplotlib
import numpy as np
import warnings

import stem_pytools.STEM_parsers as sp
from mpl_toolkits import basemap


class STEM_Domain(object):

    def __init__(self, fname_topo=None):
        """class constructor for STEM_Domain.  Populates fields: fname_topo,
        bnd_lon, bnd_lat.
        """

        self.fname_topo = fname_topo
        if self.fname_topo is None:
            self.fname_topo = os.path.join(os.getenv('SARIKA_INPUT'),
                                           'TOPO-124x124.nc')
        self.get_STEM_perimeter_latlon()
        return(None)

    def get_STEM_perimeter_latlon(self):
        """Given a STEM topo file in self.fname_topo with latitude and
        longitude coordinates of the STEM grid, populates self.bnd_lon
        and self.bnd_lat with 1-D numpy arrays containing the
        coordinates of the perimiter grid cells.
        """

        self.STEM_lon, self.STEM_lat, self.topo = (
            sp.parse_STEM_coordinates(self.fname_topo))
        self.bnd_lon = get_2d_perimeter(self.STEM_lon)
        self.bnd_lat = get_2d_perimeter(self.STEM_lat)

        return(self)

    def get_mapobj(self):
        """Return a mpl_toolkits.basemap.Basemap object for a Azimuthal
        Equidistant projection centered on the North Pole.  This
        projection accommodates in-domain tests for all latitudes
        north of 60 degrees South without encountering "wrapping"
        around the edge of the projection.
        """
        return(basemap.Basemap(projection='aeqd', lat_0=90, lon_0=0))

    def get_lat(self):
        """return the STEM latitude grid from self.topo_fname"""
        return(self.STEM_lat)

    def get_lon(self):
        """return the STEM longitude grid from self.topo_fname"""
        return(self.STEM_lon)

    def get_topo(self):
        """return the STEM topo grid from self.topo_fname"""
        return(self.topo)

    def in_domain(self, lon, lat):
        """Given numpy arrays (any dimensions) lon and lat containing
        longitude and latitude coordinates returns a boolean array of
        the same dimensions as lon and lat containing true if the
        corresponding point is inside of the STEM domain (defined as
        the area enclosed by the grid in self.fname_topo).
        """
        m = self.get_mapobj()
        x, y = m(self.bnd_lon, self.bnd_lat)
        # place projected boundary coordinates in Nx2 array
        bnd_array = np.vstack((x, y)).transpose()
        p = matplotlib.path.Path(bnd_array)

        is_in = np.zeros(lon.shape, dtype=bool)
        for i in range(lon.shape[0]):
            for j in range(lon.shape[1]):
                is_in[i, j] = p.contains_point(m(lon[i, j],
                                                 lat[i, j]))
        return(is_in)


def get_2d_perimeter(arr):
    """return a 1-D numpy array containing the perimeter values of a
    2D numpy array, starting in the top left corner and proceeding
    clockwise. This is useful for constructing a boundary path toward
    the goal of a point-in-polygon test.

    INPUTS
       arr: a 2-D numpy array

    OUTPUTS
       a 1-D numpy array containing the perimeter values of arr

    EXAMPLE
    for input arr of

    array([[0, 1, 2],
           [3, 4, 5],
           [6, 7, 8]]),

    returns:
        array([0, 1, 2, 2, 5, 8, 6, 7, 8, 6])

    """
    return(np.concatenate((arr[0, :],
                           arr[:, -1],
                           arr[-1, ::-1],
                           arr[::-1, 0])))


def calc_grid_area(xlon1, xlon2, ylat1, ylat2):
    """
    ----------------------------------------------------------------------
    DESC: calculate grid area as a function of longitude and latitude
    ----------------------------------------------------------------------
    INPUT:
    xlon1 - starting longitude
    xlon2 - ending longitude
    xlat1 - starting latitude
    xlat2 - ending latitude
    ----------------------------------------------------------------------
    OUTPUT:
    garea - grid area in m2
    ----------------------------------------------------------------------
    HISTORY:
    By: Y. F. Cheng (yafang.cheng@gmail.com)
    On: April 23, 2010

    ported from Fortran to Python by Timothy W. Hilton (thilton@ucmerced.edu)
    On: 17 Apr 2015
    ----------------------------------------------------------------------
    """
    EARTH_RADIUS = 6370000.0
    Rad_per_Deg = np.pi / 180.0

    ylat1 = fix_pole(ylat1)
    ylat2 = fix_pole(ylat2)

    cap_ht = EARTH_RADIUS * (1 - np.sin(ylat1 * Rad_per_Deg))
    cap1_area = 2 * np.pi * EARTH_RADIUS * cap_ht

    cap_ht = EARTH_RADIUS * (1 - np.sin(ylat2 * Rad_per_Deg))
    cap2_area = 2 * np.pi * EARTH_RADIUS * cap_ht

    garea = abs(cap1_area - cap2_area) * abs(xlon1 - xlon2) / 360.0

    if (np.any(abs(garea) < 1e-6)):
        warnings.warn('negative area calculated: {}'.format(
            np.argmin(garea)))

        print('{}, {}, {}'.format(xlon2, xlon1, (xlon2 - xlon1)))
        print('{}, {}, {}'.format(ylat2,
                                  ylat1,
                                  (np.sin(ylat2 * Rad_per_Deg) -
                                   np.sin(ylat1 * Rad_per_Deg))))

    return(garea)


def fix_pole(lat):
    """
    make sure that latitudes are within [-90, 90]
    """
    lat = np.maximum(lat, -90)
    lat = np.minimum(lat, 90)
    return(lat)


def demo():
    """Demonstrates using the domain module to find the points on a 1
    degree by 1 degree Northern Hemisphere grid that are inside the
    area specifed by a STEM topo file.  Plots a map showing the STEM
    domain as a blue box and the 1 degree by 1 degree grid NH grid as
    open circles. Every grid point inside the STEM domain is marked
    with a red x.
    """
    sib_lon, sib_lat = np.meshgrid(np.arange(-180, 180, 1),
                                   np.arange(0, 90, 1))
    domain = STEM_Domain()
    in_STEM_domain = domain.in_domain(sib_lon, sib_lat)
    mymap = domain.get_mapobj()
    mymap.drawcoastlines()
    mymap.plot(domain.bnd_lon, domain.bnd_lat, latlon=True)
    mymap.scatter(sib_lon,
                  sib_lat,
                  latlon=True,
                  c='gray',
                  facecolors='none',
                  marker='o')
    mymap.scatter(sib_lon[in_STEM_domain],
                  sib_lat[in_STEM_domain],
                  latlon=True,
                  c='red',
                  marker='x')
