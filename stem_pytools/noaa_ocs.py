"""
noaa_ocs
========

Module containing tools to parse, process, and visualize NOAA airborne
atmospheric species concentration observations.

Contains the class NOAA_OCS, which implements most of the
functionality, and the following functions:

get_STEMZ_height_ASL: calculate STEM Z levels from WRF heights and
STEM grid topography.
"""

import os.path
import re
import numpy as np
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
from brewer2mpl import qualitative as brewer_qualitative
import glob

import domain
import na_map


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
                 remove_poor_quality=True):
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
        ln = "#"  # initialize to empty comment
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

        self.obs['x_stem'], self.obs['y_stem'] = domain.find_nearest_stem_xy(
            self.obs.sample_longitude.values,
            self.obs.sample_latitude.values,
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
        d = domain.STEM_Domain(topo_fname)
        d.get_STEMZ_height(wrfheight_fname)
        n_obs = self.obs.shape[0]
        self.obs['z_stem'] = np.NaN  # initialize z to NaN
        # import pdb; pdb.set_trace()
        for i in range(n_obs):
            bin_idx = np.digitize(
                x=self.obs['sample_altitude'].values[np.newaxis, i],
                bins=d.asl[...,
                           self.obs.x_stem.values[i],
                           self.obs.y_stem.values[i]])
            # "+1" because python indices are 0-based, STEM indices are 1-based
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
        # lamda function to convert np.timedelta64 to hours
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

    def get_sites_lats_lons(self):
        """calculate the mean latitude and longitude for each unique site
        code in the data set.

        RETURNS:
           data frame containing columns ['sample_latitude',
              'sample_longitude'] andindex 'sample_site_code'
        """
        agg_vars = ['sample_latitude', 'sample_longitude', 'sample_site_code']
        data_agg = self.obs[agg_vars].groupby(
            ['sample_site_code']).aggregate(np.mean)
        return(data_agg)

    def plot_obs_site_locations(self):
        """
        plot a map of N America with obsertion sites labeled. The mean
        longitude and latitude for each unique site code in the data
        set is plotted to a map of North America.

        OUTPUTS:
        matplotlib.pyplot.figure containing the map
        """

        data_agg = self.get_sites_lats_lons()

        fig = plt.figure(figsize=(12, 12))
        ax = plt.axes()
        location_map = na_map.NAMapFigure(t_str='NOAA airborne obs locations',
                                          map_axis=ax)

        col = brewer_qualitative.Dark2['max'].mpl_colors[0]
        for this_site in np.unique(data_agg.index):
            x, y = location_map.map(data_agg.loc[this_site].sample_longitude,
                                    data_agg.loc[this_site].sample_latitude)
            plt.text(x, y, this_site,
                     color=col, size=24.0,
                     horizontalalignment='left',
                     verticalalignment='bottom')
            plt.scatter(x, y, marker='x', color=col)
            plt.scatter(x, y, marker='o',
                        edgecolors=col,
                        facecolors='none')
        return(location_map)

    def calculate_OCS_daily_vert_drawdown(self,
                                          altitude_bins=[0, 2000,
                                                         3400, np.inf],
                                          topo_fname=None,
                                          wrfheight_fname=None):
        """
        calculate daily mean vertical drawdown of OCS by site.  As per
        conversation with Elliott of 30 Oct 2014 this function defines
        drawdown as mean [OCS] 0 m to 2000 m minus mean [OCS] above 4000
        m.

        INPUTS:
        data: noaa_ocs.NOAA_OCS object
        altitude bins: four-element list defining top and bottom of
            "surface" altitude bin and "high" altitude bin for
            drawdown calculation.  Units are meters above ground
            level.  Default is [0, 2000, 3400, np.inf].

        OUTPUT
        pandas.DataFrame containing daily OCS drawdown by site
        """
        if topo_fname is None:
            topo_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                      'TOPO-124x124.nc')
        if wrfheight_fname is None:
            wrfheight_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                           'wrfheight-124x124-22levs.nc')
        d = domain.STEM_Domain(fname_topo=topo_fname)
        d.get_STEMZ_height(wrfheight_fname)

        # calulate noaa observation heights above ground level (m)
        # ground level altiude above sea level
        self.obs['alt_gl'] = d.asl[0, self.obs.x_stem, self.obs.y_stem]
        # sample altitude above ground level
        self.obs['sample_altitude_agl'] = (self.obs.sample_altitude -
                                           self.obs.alt_gl)
        self.obs['alt_bin'] = pd.cut(x=self.obs.sample_altitude_agl,
                                     bins=altitude_bins)

        self.obs['date'] = [t.to_datetime().date() for t in self.obs.datet]
        agg_vars = ['date', 'sample_site_code', 'alt_bin', 'analysis_value']
        grps = self.obs[agg_vars].groupby(['alt_bin',
                                           'date',
                                           'sample_site_code'])

        # take the mean within each date-site-altitude combination
        ocs_mean = grps.aggregate(np.nanmean)
        # remove date-site-altitude combinations with no observations,
        # index by site and date only, and sort so that all levels of a
        # site-date combinations are consecutive
        ocs_mean = ocs_mean.dropna().reset_index(0).sort()

        # rearrange the data frame to put "low" (0 to 2000 m) observations
        # and "hi" (> 4000 m) in consecutive columns.  Fill in missing
        # observations with NaN.
        lo = ocs_mean[ocs_mean.alt_bin == '({}, {}]'.format(altitude_bins[0],
                                                            altitude_bins[1])]
        hi = ocs_mean[ocs_mean.alt_bin == '({}, inf]'.format(altitude_bins[2])]
        mgd = pd.merge(lo, hi,
                       how='outer',
                       left_index=True,
                       right_index=True,
                       suffixes=('_lo', '_hi'))
        # create a new data frame indexed by site & date containing
        # OCS drawdown
        dd_df = pd.DataFrame(mgd.analysis_value_hi - mgd.analysis_value_lo,
                             columns=['ocs_dd'])

        return(dd_df)

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


def get_all_NOAA_airborne_data(noaa_dir):
    """
    parse all NOAA airborne OCS observation files from a directory into
    one NOAA_OCS object.  The OCS observation files are located by
    searching for files named *ocs*.txt within the specified
    directory.

    ARGS:
    noaa_dir: full path to a directory containing NOAA OCS observation
        files

    RETURNS:
    noaa_ocs.NOAA_OCS object containing all the parsed observations

    author: Timothy W. Hilton, UC Merced <thilton@ucmerced.edu>
    """
    all_files = glob.glob(os.path.join(noaa_dir, '*ocs*.txt'))
    data_list = [NOAA_OCS.parse_file(f) for f in all_files]
    # Now combine all of the DataFrames into one large DataFrame
    data = NOAA_OCS(obs=pd.concat([this_data.obs for
                                   this_data in data_list]))

    return(data)


def get_sites_summary(noaa_dir):
    """create a pandas data frame with site code, lon and lat for all
    sites with data files in the specified directory

    ARGS:
    noaa_dir (string): full path to a directory containing NOAA OCS observation
        files

    RETURNS:
    pandas DataFrame object with columns site_code, lon, lat
    """

    all_sites_df = get_all_NOAA_airborne_data(noaa_dir)
    summary_df = all_sites_df.obs.groupby('sample_site_code').mean()
    summary_df = summary_df[['sample_longitude', 'sample_latitude']]
    summary_df = summary_df.reset_index()
    summary_df.rename(columns={k: k.replace('sample_', '')
                               for k in summary_df.columns.values},
                      inplace=True)
    return(summary_df)
