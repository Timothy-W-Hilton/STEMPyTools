"""
This module implements inverse distance weighted nearest neighbor
interpolation.

Input and output coordinates may be specified using either spherical
or Cartesian coordinates (via two different functions).  Nearest
neighbors are found efficiently using a K-dimensional tree.

This module provides the following functions:

cartesian_idw_kdtree_interp: interpolate Cartesian coordinates using a
   K-dimensional tree and inverse distance weighed interpolation.

spherical_idw_kdtree_interp: interpolates spherical coordinates using
   cartesian_idw_kdtree_interp

lon_lat_to_cartesian: calculate three-dimensional Cartesian
   coordinates from spherical coordinates

Some notes on the approach taken here: This module interpolates in
horizontal space.  It does not interpolate in vertical space, or in
time.  It is assumed that the final two dimensions of the coordinate
arrays and the data array provided specify horizontal coordinates, and
that any other dimensions (vertical space, time, etc) precede these.
I decided to set it up this way because this module is being developed
for atmospheric chemical transport modeling, and the quantities of
interest in that context (species concentrations, wind velocities,
pressures, etc.) vary on vastly different length scales horizontally
than they do vertically.  Nearest neighbor interpolation thus becomes
a more complex problem if vertical levels are involved.  For example,
consider a model grid on a relatively coarse spatial grid (1 degree by
1 degree or something).  The vertical levels could easily be such that
the 10 nearest neighbors to any grid point are the other points in the
same vertical column.  These generally aren't the points that should
be interpolated in a meaningful regridding exercise.
"""

import numpy as np
from scipy.spatial import cKDTree

def spherical_idw_kdtree_interp(lon_s, lat_s, lon_d, lat_d,
                                data,
                                R=6371007.181000,
                                n_nbr=9):
    """
    Interpolate values on a source grid to a destination grid using
    inverse distance weighted interpolation.

    The source grid and destination grid are specified in spherical
    coordinates.

    PARAMETERS
    ==========
    lon_s, lat_s; np.ndarray in 2 dimensions: source longitudes and
        latitudes respectively.  lon_s and lat_s must have identical
        shape.
    lon_d, lat_d in 2 dimensions: destination longitudes and
        latitudes, respectively.  lon_d and lat_d must have identical
        shape.
    data: np.ndarray: the values to be interpolated.  The last 2
        dimensions of data must match the last 2 dimensions of lon_s
        and lat_s.
    R: radius of the sphere.  The default of 6371007.181000 is the
        radius of the Earth in meters.
    n_nbr: number of nearest neighbors to use in the inverse distance
        weighted interpolation.

    RETURNS
    =======
    data_idw: Numpy ndarray; the values from data interpolated onto
        the destination points using inverse distance weighted
        interpolation.  Has same number of dimensions as data, with
        the last 2 dimensions matching the dimensions of lon_d and
        lat_d.

    Notes: I don't think the radius of the sphere really matters --TWH
    """
    # if lon_s.ndim < data.ndim:
    #     lon_s = np.resize(lon_s, data.shape)
    #     lat_s = np.resize(lat_s, data.shape)
    xs, ys, zs = lon_lat_to_cartesian(lon_s,
                                      lat_s,
                                      R=R)
    xd, yd, zd = lon_lat_to_cartesian(lon_d,
                                      lat_d,
                                      R=R)
    data_idw = cartesian_idw_kdtree_interp(
        xs, ys, zs,
        xd, yd, zd,
        data,
        n_nbr=n_nbr)

    return(data_idw)

def cartesian_idw_kdtree_interp(xs, ys, zs, xd, yd, zd, data, n_nbr=9):
    """
    Interpolate values from source coordinates to destination
    coordinates using inverse distance weighted interpolation.

    xs, ys, zs; np.ndarray in 2 dimensions: source X, Y and Z
        coordinates, respectively. zs, ys, and zs must have identical
        shapes.
    xd, yd, zd; np.ndarray in 2 dimensions: destination X, Y and Z
        coordinates, respectively. zd, yd, and zd must have identical
        shapes.
    data: np.ndarray: the values to be interpolated.  The last 2
        dimensions of data must match the last 2 dimensions of xs, ys,
        and zs.
    n_nbr: number of nearest neighbors to use in the inverse distance
        weighted interpolation.

    RETURNS
    =======
    data_idw: Numpy ndarray; the values from data interpolated onto
        the destination coordinates using inverse distance weighted
        interpolation.  Has same number of dimensions as data, with
        the last 2 dimensions matching the dimensions of lon_d and
        lat_d.
    """
    # construct a K-dimensional tree of source coordinates and search
    # it for nearest neighbors of destination coordinates
    dists, inds = find_nn_3D_kdtree(xs, ys, zs, xd, yd, zd, n_nbr)

    # if n_nbr is 1, dists and inds are returned as a vector [shape is
    # (M,)], whereas for n_nbr > 1 dists and inds are returned as M by
    # n_nbr arrays [shape is (M,n_nbr)].  The np.sum axis keyword
    # below depends on that second dimension existing.  If it does not
    # exist, promote the M-element vector to an M by 1 array.
    if dists.ndim == 1:
        dists = dists[..., np.newaxis]
        inds = inds[..., np.newaxis]
    wgts = 1.0 / dists**2
    data_idw = np.zeros((data.shape[0], xd.shape[-2], xd.shape[-1]))
    print 'starting loop'
    for i in range(data.shape[0]):
        data_idw[i, ...] = (
            np.sum(wgts * data[i, ...].flatten()[inds], axis=1) /
            np.sum(wgts, axis=1)).reshape(xd.shape)

    #the interpolated data should have the shape of the input data in
    #all but the last two dimensions, and the shape of the destination
    #coordinates in the last two
    return(data_idw)

def find_nn_3D_kdtree(xs, ys, zs, xd, yd, zd, n_nbr):
    """
    Find the n_nbr nearest neighbors to one set of coordinates from
    within a second set of coordinates.

    A K-dimensional tree is constructed from a set of source coordinates
    and searched for the n_nbr nearest neighbors to the points in the
    destination coordinates.  All coordinates are Cartesian.

    PARAMETERS
    ==========
    xs, ys, zs; np.ndarray in 2 dimensions: source X, Y and Z
        coordinates, respectively. zs, ys, and zs must have identical
        shapes.
    xd, yd, zd; np.ndarray in 2 dimensions: destination X, Y and Z
        coordinates, respectively. zd, yd, and zd must have identical
        shapes.
    data: np.ndarray: the values to be interpolated.  The last 2
        dimensions of data must match the last 2 dimensions of xs, ys,
        and zs.
    n_nbr: scalar integer; number of nearest neighbors to find

    RETURNS
    dists: np.ndarray, xd.size by n_nbr; the distances from each
       destination point to its n_nbr nearest neighbors specified by
       inds.
    inds: np.ndarray, xd.size by n_nbr; indices the n_nbr nearest neighbors
       to each destination point.  The indices are flat indices.

    """
    tree = cKDTree(zip(xs.flatten(), ys.flatten(), zs.flatten()))
    dists, inds = tree.query(zip(xd.flatten(),
                                 yd.flatten(),
                                 zd.flatten()),
                             n_nbr)
    return dists, inds

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    Convert spherical coordinates to Cartesian coordinates on a sphere
    with radius R.

    PARAMETERS
    ==========
    lon, lat; np.ndarray: longitudes and latitudes of the points.  lon
       and lat must have the same shape.
    R: scalar; the radius of the sphere

    RETURNS
    x, y, z: np.ndarray of same shape as lon, lat: The Cartesian
        coordinates of the points.

    Documented by Timothy W. Hilton. Written and posted at earthpy.org
    by Oleksandr Huziy.
    http://earthpy.org/interpolation_between_grids_with_ckdtree.html
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z    
