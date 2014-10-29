"""
noaa_ocs
========

Module containing tools to parse, process, and visualize NOAA airborne
atmospheric species concentration observations.

Contains the class NOAA_OCS, which implements most of the
functionality, and the following functions:

get_STEMZ_height_ASL: calculate STEM Z levels from WRF heights and
STEM grid topography.

find_nearest_stem_xy: find the STEM grid cell nearest to a specified
latitude and longitude.

lon_lat_to_cartesian: calculate cartesian X, Y, Z coordinates from
spherical coordinates.
"""

import netCDF4
import re
import math, os
import numpy as np
from datetime import datetime
import pandas as pd
import csv
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
from scipy.spatial import cKDTree
from brewer2mpl import qualitative as brewer_qualitative

import na_map
import STEM_parsers


class NOAA_OCS(object):
    """
    Container class to work with NOAA airborne observations
    
    Class attributes:
       obs: pandas DataFrame containing NOAA observations
       obs_color: matplotlib color to use for plotting observations.  
       grid_color: matplotlib color to use for plotting STEM grid cell points
    """
    def __init__(self,
                 obs=None,
                 remove_poor_quality = True):
        """
        Class constructor; builds object from an existing Pandas
        DataFrame of observations.

        May be called directly with an existing DataFrame, but
        typically it will be most useful to use the class method
        parse_file to create a NOAA_OCS object from a NOAA airborne
        obseration CSV file.
        
        PARAMETERS
        ----------
        obs: Pandas DataFrame containing NOAA airborne observations.
            There are no limit on the number of observations or the
            number of variables, but the following variables are
            required: sample_longitude, sample_latitude, sample_year,
            sample_month, sample_day, sample_hour, sample_minute,
            sample_seconds.
        """
        self.obs = obs
        self.obs_color = brewer_qualitative.Dark2['max'].mpl_colors[0]
        self.grid_color = brewer_qualitative.Dark2['max'].mpl_colors[1]
        if remove_poor_quality:
            self = self.remove_poor_quality_obs()

    @classmethod
    def parse_file(cls, fname):
        """
        Construct a NOAA_OCS object from a data file.

        Parse a NOAA airborne measurement whitespace-delimited ASCII
        file and create a NOAA_OCS object containing the observations
        in its obs field. The last comment line in the file header
        should contain the column names in the format "data_fields:
        NAME_1 NAME_2 NAME_3 ... NAME_N".  There are no restrictions
        on the number of variables in the file, but the following
        variables are required: sample_longitude, sample_latitude,
        sample_year, sample_month, sample_day, sample_hour,
        sample_minute, sample_seconds.  The timestamp is converted to
        a datetime.datetime object and placed in the obs DataFrame as
        obs.datet.
        PARAMETERS
        ----------
        fname: full path the NOAA file to be parsed

        returns: object of class noaa_ocs.NOAA_OCS.
        
        """

        # the NOAA files contain varying number of header lines.  It
        # looks like the column names are always the last commented
        # header line.  Now get the column names and the number of 
       # header lines.
        f = open(fname, 'r')
        ln = "#"  #initialize to empty comment
        hdr_count = 0
        while(ln[0] is '#'):
            last_ln = ln
            ln = f.readline()
            hdr_count += 1
        names = re.sub('.*data_fields:\s+', '', last_ln).rstrip().split()
        # parse the data to a DataFrame
        df = pd.read_csv(fname,
                         comment='#',
                         delim_whitespace=True,
                         skiprows=hdr_count-1,
                         names=names)
        # convert the times to datetime objects
        df['datet'] = df.apply(
            lambda row: datetime(row['sample_year'],
                                 row['sample_month'],
                                 row['sample_day'],
                                 row['sample_hour'],
                                 row['sample_minute'],
                                 row['sample_seconds']),
            axis=1)

        return(cls(obs=df))

    def remove_poor_quality_obs(self):
        """
        Remove observations with a "do not use" QC flag.  "Do not use"
        is currently defined as any character other than '.' in the
        first of the three QC columns.

        From the NOAA OCS data readme file:
        -----
        NOAA ESRL uses a 3-column quality control flag where each
        column is defined as follows:

        column 1    REJECTION flag.  An alphanumeric other
                    than a period (.) in the FIRST column indicates
                    a sample with obvious problems during collection
                   or analysis.  This measurement should not be interpreted.
        -----

        The second and third column are flags set for the needs of
        specific PIs or projects.  For now, at least, the noaa_ocs
        module ignores these.
        """
        quality_good = self.obs.apply(
            lambda row: row.analysis_flag[0] == '.',
            axis=1)
        self.obs = self.obs[quality_good]

    def get_stem_xy(self, stem_lon, stem_lat):
        """
        Get STEM grid (x, y) values for observations.

        Uses a nearest neighbor search to find the STEM grid cell
        center nearest the observation location.  Places the STEM x
        and y grid indices in the variables obs.x_stem and obs.y_stem,
        respectively.  Places the STEM longitude and latitude in
        variables obs.lon_stem and obs.lat_stem, respectively.
        PARAMETERS
        ----------
        stem_lon: longitudes of STEM grid cell centers
        stem_lat: latitudes of STEM grid cell centers
        """
        self.obs['x_stem'], self.obs['y_stem'] = find_nearest_stem_xy(
            self.obs.sample_longitude,
            self.obs.sample_latitude,
            stem_lon,
            stem_lat)
        self.obs['lon_stem'] = stem_lon[self.obs.x_stem, self.obs.y_stem]
        self.obs['lat_stem'] = stem_lat[self.obs.x_stem, self.obs.y_stem]

    def get_stem_z(self,
                   topo_fname='./TOPO-124x124.nc',
                   wrfheight_fname='./wrfheight-124x124-22levs.nc'):
        """
        Get STEM grid z values for observations.

        Determine the STEM z level of NOAA airborne observations using
        each observation's STEM grid x, STEM grid y, and altitude
        above sea level, and STEM grid sigma levels.  Places the STEM
        Z level in the variable obs.z_stem.
        PARAMETERS
        ----------
        topo_fname: full path the to the IO/API file defining the
           latitude, longitude, and topopgraphy of the STEM grid
        wrfheight_fname: full path the to the IO/API file specifying
           the WRF heights on the STEM grid.
        """
        asl = get_STEMZ_height_ASL(topo_fname, wrfheight_fname)
        n_obs = self.obs.shape[0]
        self.obs['z_stem'] = np.NaN  #initialize z to NaN
        #import pdb; pdb.set_trace()
        for i in range(n_obs):
            bin_idx = np.digitize(
                x=self.obs['sample_altitude'].values[np.newaxis, i],
                bins=asl[...,
                         self.obs.x_stem.values[i],
                         self.obs.y_stem.values[i]])
            #"+1" because python indices are 0-based, STEM indices are 1-based
            self.obs['z_stem'].values[i] = bin_idx + 1 

    def get_stem_t(self, stem_t0):
        """
        get STEM timestep for airborne observations.

        Calculate the stem timestep (in hours from STEM start time)
        for each observation.  The STEM timestep is placed in the obs
        DataFrame variable stemt.
        PARAMETERS
        ----------
        stem_t0: datetime.datetime; the first timestep of the stem
           simulation.
        """
        #lamda function to convert np.timedelta64 to hours
        td2hrs = lambda x: x / np.timedelta64(1, 'h')
        self.obs['stemt'] = (self.obs.datet - stem_t0).apply(td2hrs)

    def plot_stem_obs_matches(self):
        """
        Plot observations and their assigned grid cells on a map.

        Plot observations (^) and nearest STEM grid cells (o) on a map
        of North America, with a line connecting each observation to
        its assigned grid cell.
        """
        lats = np.array((self.obs.sample_latitude.values,
                         self.obs.lat_stem.values))
        lons = np.array((self.obs.sample_longitude.values,
                         self.obs.lon_stem.values))
        lat_min = lats.min()
        lat_max = lats.max()
        lon_min = lons.min()
        lon_max = lons.max()
        lat_pad = (lat_max - lat_min) * 0.1
        lon_pad = (lon_max - lon_min) * 0.1

        # lon, lat for SW corner of map
        SW_crnr = (max(lon_min - lon_pad, -180),
                   max(lat_min - lat_pad, -90)) 
        # lon, lat for NE corner of map
        NE_crnr = (min(lon_max + lon_pad, 720),
                   min(lat_max + lat_pad, 90))

        m = na_map.NAMapFigure(cb_axis=None).map
        # m = basemap.Basemap(llcrnrlon=SW_crnr[0], llcrnrlat=SW_crnr[1],
        #                     urcrnrlon=NE_crnr[0], urcrnrlat=NE_crnr[1])
        # m.drawmeridians(np.linspace(lon_max, lon_min, num=3),
        # 			labels=[0,0,0,1])
        # m.drawparallels(np.linspace(lat_max, lat_min, num=3),
        # 			labels=[1,0,0,0])

        obs_mrks = m.scatter(self.obs.sample_longitude.values,
                             self.obs.sample_latitude.values,
                             c=self.obs_color,
                             label='observation',
                             zorder=2,
                             latlon=True)
        grid_mrks = m.scatter(self.obs['lon_stem'].values,
                              self.obs['lat_stem'].values,
                              color=self.grid_color,
                              marker='^',
                              label='STEM grid cell',
                              zorder=3,
                              latlon=True)
        n_obs = self.obs.shape[0]
        for i in range(n_obs):
            ln = m.plot((self.obs['lon_stem'].values[i],
                         self.obs['sample_longitude'].values[i]),
                        (self.obs['lat_stem'].values[i],
                         self.obs['sample_latitude'].values[i]),
                        color='black',
                        latlon=True,
                        zorder=1)
        plt.legend(numpoints=1,
                   scatterpoints=1,
                   loc='best')
        return(m)

    def plot_locations(self):
        """
        Draw observation locations on a map.

        Draw a scatter plot on a map of North America of the location
        of each observation.
        """

        m = na_map.NAMapFigure(cb_axis=None).map
        mrks = m.scatter(
            self.obs.sample_longitude.values,
            self.obs.sample_latitude.values,
            c=self.obs_color,
            label='observation',
            latlon=True)
        return(m, mrks)

# def init_NA_map_cyl():
#     """
#     Draw a map of North America using a cylindrical projection.
#     """
#     plt.figure(figsize=(8,8))
#     #initialize a map of North America
#     SW_crnr = (-170, 10) # lon, lat for SW corner of map
#     NE_crnr = (-50, 80) # lon, lat for NE corner of map
#     m = basemap.Basemap(llcrnrlon=SW_crnr[0], llcrnrlat=SW_crnr[1],
#                         urcrnrlon=NE_crnr[0], urcrnrlat=NE_crnr[1])
#     m.drawmeridians(np.linspace(SW_crnr[0], NE_crnr[0], num=3),
#                     labels=[1,0,0,0])
#     m.drawparallels(np.linspace(SW_crnr[1], NE_crnr[1], num=3),
#                     labels=[0,0,0,1])
#     m.drawcoastlines()
#     return(m)

def get_STEMZ_height_ASL(topo_fname='./TOPO-124x124.nc',
                         wrfheight_fname='./wrfheight-124x124-22levs.nc'):
    """
    Get STEM Z level boundaries in meters above sea level.

    Read the wrf height and topo files to get the STEM vertical cell
    boundaries in meters above sea level.  How model levels should be
    interpreted is described in (1) and (2):
    (1) http://www.ecmwf.int/research/ifsdocs/DYNAMICS/Chap2_Discretization4.html#961180
    (2) Simmons, A.J., and D.M. Burridge, 1981: An energy and angular
        momentum conserving vertical finite difference scheme and hybrid
        vertical coordinates. Mon. Weather Rev., 109, 758-766.
    PARAMETERS
    ----------
    topo_fname: full path to the IO/API file containing lat, lon, and
       height above sea level for STEM grid cells
    wrfheight_fname: full path to the IO/API file containing z levels
       for stem grid cells
    """

    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(topo_fname)

    nc = netCDF4.Dataset(wrfheight_fname, 'r', format='NETCDF4')
    #numpy array shape(22, 124, 124): zlev, lat, lon.  squeeze removes
    #the time dimension
    agl = nc.variables['AGL'][:].squeeze()
    nc.close()

    asl = agl + topo #[np.newaxis, ...].shape
    return(asl)

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
    a 2-element tuple of X and Y indices, one index per observation.
    """
    # convert spherical lon, lat coordinates to cartesian coords. Note
    # that these x,y,z are 3-dimensional cartesian coordinates of
    # positions on a sphere, and are thus different from the x,y,z
    # *indices* of the STEM grid.
    x, y, z = lon_lat_to_cartesian(lon, lat)
    xs, ys, zs = lon_lat_to_cartesian(lon_stem, lat_stem)

    #use a K-dimensional tree to find the nearest neighbor to (x,y,z)
    #from the points within (xs, ys, zs).  A KD tree is a data
    #structure that allows efficient queries of K-dimensional space (K
    #here is 3).
    tree = cKDTree(np.dstack((xs.flatten(),
                              ys.flatten(),
                              zs.flatten())).squeeze())
    d, inds = tree.query(
        np.dstack((x.values,
                   y.values,
                   z.values)).squeeze(), k = 1)

    return(np.unravel_index(inds, lon_stem.shape))

def lon_lat_to_cartesian(lon, lat, R = 1):
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

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z


