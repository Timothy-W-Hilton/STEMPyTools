"""a set of tools useful for parsing and returning descriptions of the
STEM domain (grid cell longitude and latitude coordinates, grid cell
vertical coordinates), and also for some domain-oriented calculations
(is a point inside or outside of the domain, etc.)
"""

import os
import os.path
import matplotlib
import numpy as np
import warnings
import netCDF4
import stem_pytools.STEM_parsers as sp
from scipy.spatial import cKDTree


class STEM_Domain(object):
    """class to contain STEM domain attributes.
    """
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

    def get_lat(self):
        """return the STEM latitude grid from self.fname_topo"""
        return(self.STEM_lat)

    def get_lon(self):
        """return the STEM longitude grid from self.fname_topo"""
        return(self.STEM_lon)

    def get_topo(self):
        """return the STEM topo grid from self.fname_topo"""
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

    def get_STEMZ_height(self,
                         wrfheight_fname=None):
        """
        Get STEM Z level boundaries in meters above sea level.

        Read the wrf height and topo files to get the STEM vertical cell
        boundaries in meters above sea level.  How model levels should be
        interpreted is described in (1) and (2):
        (1)
        http://www.ecmwf.int/research/ifsdocs/DYNAMICS/Chap2_Discretization4.html#961180
        (2) Simmons, A.J., and D.M. Burridge, 1981: An energy and angular
            momentum conserving vertical finite difference scheme and hybrid
            vertical coordinates. Mon. Weather Rev., 109, 758-766.
        PARAMETERS
        ----------
        topo_fname: full path to the IO/API file containing lat, lon, and
           height above sea level for STEM grid cells
        wrfheight_fname: full path to the IO/API file containing z levels
           for stem grid cells

        OUTPUT
        two-element tuple of numpy ndarrays: (agl, asl).  agl contains the
        stem cell heights above ground level, asl contains the stem cell
        heights above sea level.
        """
        if wrfheight_fname is None:
            wrfheight_fname=os.path.join(os.getenv('SARIKA_INPUT'),
                                         'wrfheight-124x124-22levs.nc')

        # numpy array shape(22, 124, 124): zlev, lat, lon.  squeeze removes
        # the time dimension
        nc = netCDF4.Dataset(wrfheight_fname, 'r', format='NETCDF4')
        self.agl = nc.variables['AGL'][:].squeeze()
        nc.close()
        # parse_STEM_var currently won't handle time-independent I/O
        # API files
        # agl_var = sp.parse_STEM_var(wrfheight_fname,
        #                             varname='AGL')
        # self.agl = agl_var['data'].squeeze()
        self.asl = self.agl + self.get_topo()


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


def get_stem_z_from_altitude(alts,
                             stem_x=None,
                             stem_y=None,
                             lons=None,
                             lats=None,
                             topo_fname=None,
                             wrfheight_fname=None):
    """
    Get STEM grid z values for specified altitudes

    Determine the STEM z level of specified altutides using each
    observation's STEM grid x, STEM grid y, and altitude above sea
    level, and STEM grid sigma levels.  Either (stem_x and stem_y) or
    (lons and lats) must be specified.  If stem_x, stem_y, lons and
    lats are all specified stem_x and stem_y will be used and lons and
    lats ignored.  alts, stem_x, stem_y, lons and lats must all
    contain the same number of elements.

    ARGS:
    alts (array-like): alitude values to place in STEM Z cells (meters)
    stem_x (array-like): STEM domain x coordinates
    stem_y (array-like): STEM domain y coordinates
    lons: (array-like): longitude coordinates of alts
    lats: (array-like): latitude coordinates of alts
    topo_fname: full path the to the IO/API file defining the
       latitude, longitude, and topopgraphy of the STEM grid.
       Defaults to $SARIKA_INPUT/TOPO-124x124.nc
    wrfheight_fname: full path the to the IO/API file specifying the
       WRF heights on the STEM grid.  Defaults to
       $SARIKA_INPUT/wrfheight-124x124-22levs.nc

    RETURNS:
    array the same shape as alts containing the STEM Z levels
    """

    d = STEM_Domain(topo_fname)

    if stem_x is None or stem_y is None:
        if lons is None or lats is None:
            raise TypeError(("you must specify both (stem_x and stem_y) "
                             "OR both (lons and lats)"))
        stem_x, stem_y = find_nearest_stem_xy(lons, lats,
                                              d.get_lons(), d.get_lats())

    d.get_STEMZ_height(wrfheight_fname)
    shape_alts = alts.shape
    alts = alts.flatten()

    z_stem = np.zeros(alts.shape, dtype=int)
    z_stem[:] = -999
    for i in range(alts.size):
        bin_idx = np.digitize(
            x=alts,
            bins=d.asl[..., stem_x[i], stem_y[i]])
        # "+1" because python indices are 0-based, STEM indices are 1-based
        z_stem[i] = bin_idx[i] + 1
    z_stem = z_stem.reshape(shape_alts)
    return(z_stem)


def lon_lat_to_cartesian(lon, lat, R=1):
    """
    Convert spherical coordinates to three-dimensional Cartesian
    coordinates.

    calculates three dimensional cartesian coordinates (x, y, z) for
    specified longitude, latitude coordinates on a sphere with radius
    R.  Written and posted at earthpy.org by Oleksandr Huziy.
    http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    accessed 19 Mar 2014 by Timothy W. Hilton.

    PARAMETERS
    ==========
    lon; np.ndarray: longitude values
    lat; np.ndarray: latitude values
    R: scalar; radius of the sphere.

    RETURNS
    =======
    three element tuple containing X, Y, and Z, one element per
    lon,lat pair.
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x, y, z


def find_nearest_stem_xy(lon, lat, lon_stem, lat_stem):
    """
    find nearest neighbors for a set of lon, lat points from a second
    set of lon, lat points.

    Given a set of arbitrary (lon, lat) positions, find the horizontal
    (x, y) STEM grid indices of the nearest STEM grid cell center to
    each position.
    PARAMETERS
    ----------
    lon, lat: ndarray; of arbitrary longitudes and latitudes.  Must
       contain the same number of elements.
    lon_stem, lat_stem: ndarrays; longitudes and latitudes of STEM
       grid cell centers. Must contain the same number of elements.

    RETURNS:
    an N-element tuple of X and Y indices, one index per observation.
        The indices are the closest point in (lon, lat) to each point
        in (lon_stem, lat_stem).  N is therefore equal to the number
        of dimensions in lon_stem and lat_stem.
    """
    # convert spherical lon, lat coordinates to cartesian coords. Note
    # that these x,y,z are 3-dimensional cartesian coordinates of
    # positions on a sphere, and are thus different from the x,y,z
    # *indices* of the STEM grid.
    x, y, z = lon_lat_to_cartesian(lon, lat)
    xs, ys, zs = lon_lat_to_cartesian(lon_stem, lat_stem)

    # use a K-dimensional tree to find the nearest neighbor to (x,y,z)
    # from the points within (xs, ys, zs).  A KD tree is a data
    # structure that allows efficient queries of K-dimensional space (K
    # here is 3).
    tree = cKDTree(np.dstack((xs.flatten(),
                              ys.flatten(),
                              zs.flatten())).squeeze())

    d, inds = tree.query(
        np.dstack((x, y, z)).squeeze(), k=1)

    return(np.unravel_index(inds, lon_stem.shape))


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
