"""Provides functionality useful for postporcessing STEM AQOUT files

class aqout_container() provides most of the functionality.  Several
helper functions are also included.
"""

import sys
import numpy as np
import os.path
import os
import pandas as pd
import netCDF4
import warnings
import cPickle
from datetime import datetime, timedelta
import itertools
import subprocess

from stem_pytools import STEM_parsers as sp
from stem_pytools.check_paths import check_path_with_msg
from stem_pytools import NERSC_data_paths as ndp
from stem_pytools import calc_drawdown


class aqout_container(object):
    """class to contain and optionally combine data from one or more
    STEM AQOUT files.

    CONTAINS FIELDS:
    aqout_paths: list or tuple; full paths to the object's AQOUT file(s)
    desc: character string; text describing the COS scenario described
       in the object's data
    key: character string; short string for use a key in a dict of
       aqout scenarios

    The class also provides a str method to produce a nicely formatted
    printout of the paths it contains.

    """
    def __init__(self,
                 aqout_paths=[],
                 desc='',
                 key=''):
        """class constructor"""
        # if the aqout_paths argument contains a single string,
        # convert to a one-element tuple.  Doing that here rather than
        # requiring that the aqout_paths be a tuple already allows for
        # creating arguments for the constructor by, for example,
        # calling itertools.product on two lists of strings.
        if type(aqout_paths) is str:
            aqout_paths = (aqout_paths,)
        self.desc = desc
        self.key = key
        self.aqout_paths = aqout_paths
        self.data = []  # list to hold data arrays
        self.t = []  # list to hold time stamp arrays
        self.cos_total = None  # field to hold cos total

    def __str__(self):
        """ formatted printing of the paths in an aqout_combiner object. """
        return('\n'.join(self.aqout_paths))

    def all_paths_exist(self):
        """Checks whether all files in object's aqout_paths exist.  If so
        returns True.  If not, writes the path(s) that do/does not
        exist to stdout and return False. """
        all_exist = all([check_path_with_msg(p) for p in self.aqout_paths])
        return(all_exist)

    def parse(self, t0, t1, verbose=False):
        """parse the concentrations in the object's AQOUT file(s) for the
        period beginning with t0 and ending with t1.  The
        concentration are placed in a list of numpy ndarray objects;
        that list is placed in the calling object's 'data' field.  The
        corresponding timestamps are also placed in a list.  That list
        is then placed in the calling object's 't' field.

        Assumes that the timestep of all AQOUT files is one hour.
        Results are undefined if this requirement is not met.

        INPUTS
        t0, t1: datetime.datetime objects specifying the starting and
            ending timestamps to combine

        """

        one_hour = 10000  # Models-3 I/O API format for 1 hour
                          # (integer, HHMMSS)

        for p in self.aqout_paths:
            nc = netCDF4.Dataset(p)
            tstep = nc.TSTEP
            nc.close()
            if tstep != one_hour:
                raise ValueError('timestep of {} is not one hour!'
                                 '(it is {})'.format(os.path.basename(p),
                                                     tstep))
            if verbose:
                sys.stdout.write('parsing {}\n'.format(os.path.basename(p)))
                sys.stdout.flush()
            this_cos = sp.parse_STEM_var(p,
                                         varname='CO2_TRACER1',
                                         t0=t0,
                                         t1=t1)
            self.data.append(this_cos['data'])
            this_t = pd.DatetimeIndex(this_cos['t'], freq='1H')
            self.t.append(this_t)

    def sum(self):
        """Add AQOUT concentration fields together.  For now uses the
        (simple-as-can-be) approach of adding together the parsed
        concentration fields element-wise.  I've implemented the
        adding as its own method so that the infrastructure is there
        in the future to handle more complicated cases (such as AQOUT
        files containing different timesteps).  """
        if len(self.data) == 0:
            raise ValueError('object {} contains no '
                             'un-summed AQOUT data'.format(self.key))

        assert all([self.t[0].equals(this_t) for this_t in self.t[1:]])
        self.t = self.t[0]

        self.cos_total = reduce(np.add, self.data)
        del self.data
        self.data = []

    def calc_stats(self):
        """

        """
        if self.cos_total is None:
            raise ValueError('object {} contains no '
                             'summed AQOUT data'.format(self.key))

        self.t_stats, self.cos_mean = daily_window_stats(self.t,
                                                         self.cos_total,
                                                         is_midday,
                                                         np.mean)
        self.t_stats, self.cos_std = daily_window_stats(self.t,
                                                        self.cos_total,
                                                        is_midday,
                                                        np.std)

    def calc_drawdown(self,
                      topo_fname=None, wrfheight_fname=None,
                      lo_height_agl=2000, hi_height_agl=4000):
        """calculate mean COS drawdown from cos_total.  Arguments are
        passed directly to calc_drawdown.calc_JA_midday_drawdown.
        """

        dd = calc_drawdown.calc_STEM_COS_drawdown(self.cos_total,
                                                  topo_fname=None,
                                                  wrfheight_fname=None)
        return dd

    def stats_to_netcdf(self, fname):

        """Write daily COS mean and standard deviation to a netCDF file.
        This is meant to replace the cPickle saving.  My thinking here
        is that these STEM results are useful, hopefully in the
        slightly longer term, and a more portable,
        platform-independent storage and transfer format is useful.

        Note: netCDF file is opened in 'write' mode -- this means an
        existing file named fname will be deleted."""

        # TO DO: check here that mean, std dev have in fact been calculated

        NC_FLOAT = 'f'  # netCDF4 specifier for NC_FLOAT datatype
        NC_INT = 'i4'  # netCDF4 specifier for NC_DOUBLE datatype

        if os.path.exists(fname):
            warnings.warn(("{} already exists. Please delete "
                           "or rename before proceeding. "
                           "Exiting; file not written".format(fname)))
            return
        nc = netCDF4.Dataset(fname, 'w')

        nc.createDimension('T', self.cos_mean.shape[0])
        nc.createDimension('LAY', self.cos_mean.shape[1])
        nc.createDimension('ROW', self.cos_mean.shape[2])
        nc.createDimension('COL', self.cos_mean.shape[3])

        nc.createVariable(varname='cos_mean',
                          datatype=NC_FLOAT,
                          dimensions=(('T', 'LAY', 'ROW', 'COL')))
        nc.createVariable(varname='cos_std',
                          datatype=NC_FLOAT,
                          dimensions=(('T', 'LAY', 'ROW', 'COL')))
        nc.createVariable(varname='time',
                          datatype=NC_INT,
                          dimensions=(('T')))

        # describe units, etc.
        nc.variables['time'].description = 'time stamp'
        nc.variables['time'].units = 'YYYYMMDDhhhmmmss'
        var = nc.variables['cos_mean']
        var.description = 'mid-day COS concentration mean'
        var.units = 'molecules cm-3'
        var = nc.variables['cos_std']
        var.description = ('mid-day COS concentration'
                           'standard deviation')
        var.units = 'molecules cm-3'
        nc.key = self.key
        nc.description = (self.desc +
                          ('July and August North '
                           'American mid-day [COS]'
                           'mean and standard deviation.  Mid-day is '
                           'defined as between 15:00 and 23:00 UTC, '
                           'which is 10:00 EST to 15:00 PST.'))

        nc.variables['cos_mean'][...] = self.cos_mean
        nc.variables['cos_std'][...] = self.cos_std
        nc.variables['time'][...] = map(
            lambda x: np.int(x.strftime('%Y%m%d%H%M%S')), self.t_stats)

        nc.close()


    def align_tstamps(self):
        """This approach will work to combine aqout files at timestamps
        that are present in both.  But it won't work to fill in timestamps
        that are only present in one.  That's OK for now because all our
        AQOUT files are hourly.  I wonder if I should plan ahead now for
        AQOUT files with a different output time step.

        In the case of different timesteps I think combining AQOUT files
        via calls to the Models-3 I/O API routines (currstep, etc) would
        be a better approach.  I wonder how difficult it would be to
        access that fortran via python...  Would certainly be quicker to
        just do the whole thing in Fortran... """

        warnings.warn('this function is a placeholder/work in progress')
        dfs = [pd.DataFrame(data=np.arange(len(t)),
                            index=t)
               for t in self.t]
        # # this works as intended
        # dfs[0].[1:8].merge(dfs[1].iloc[5:],
        #                    how='inner',
        #                    left_index=True, right_index=True)
        dfs[0].merge(dfs[1], how='inner',
                     left_index=True, right_index=True)
        return(dfs)


def is_midday(this_t):
    """
    Return true if this_t is "mid-day" in North America.  "Mid-day" is
    defined as between 15:00 and 23:00 UTC, which is 10:00 EST to
    15:00 PST.

    INPUTS
    this_t: a datetime.datetime object

    OUTPUTS
    True or False
    """
    return ((this_t.hour >= 15) and (this_t.hour <= 23))


def daily_window_stats(t, data, f_wndw, f_stat):
    """
    calculate a statistic within a daily time window.

    INPUT PARAMETERS
    data: np.ndarray of shape [t,z,x,y] with axes time, vertical
       coordinate, horizontal coordinate 1, horizontal coordinate 2.
    t: pandas DatetimeIndex of length t
    f_wndw: a function that accepts a datetime.datetime object and
       returns true if the object is inside the desired daily window and
       false if outside.
    f_stat: function to calculate the statistic.  For example np.mean
       or np.std.  Must accept "axis=0" as an argument.

    OUTPUT PARAMETERS
    data_out: np.ndarray of shape [t1, z, x, y] with the calculated
       statistic
    t_out: np.ndarray of datetime.date objects with the date of the
       data in data_out
    """
    idx_wndw = np.where(t.to_series().apply(f_wndw))[0]
    t = t[idx_wndw]
    data = data[idx_wndw, ...]

    t_date = t.to_series().apply(datetime.date)
    idx_split = np.where(t_date.diff() > timedelta(0))[0]
    lst = np.split(data, idx_split)
    data_lst = list(itertools.imap(lambda x: f_stat(x, axis=0), lst))

    data_out = np.array(data_lst)
    t_out = t_date.unique()

    return(t_out, data_out)


def assemble_data(model_runs=None, pickle_fname=None):
    """calculate and save to a cPickle file:
    (1) daily mid-day [COS] mean
    (2) daily mid-day [COS] standard deviation
    (3) timestamps associated with (1) and (2)

    the data are placed in a dict and saved to a user-specified
    cpickle file.

    INPUT PARAMETERS
    model_runs: dict of STEMRun objects.  If unspecified the default
        is the output of stem_pytools.NERSC_data_paths.get_runs()
    pickle_fname: full path of the cpickle file to create.  If
        unspecified the default is
        /home/thilton/Data/STEM/aq_out_data.cpickle

    OUTPUT PARAMETERS:
    all_data: dict containing (1), (2), and (3) above.

    """

    warnings.warn('assemble_data is deprecated and will be removed in'
                  ' the future.  Please use '
                  'aqout_postprocess.aqout_container instead')

    if model_runs is None:
        model_runs = ndp.get_runs()

    t = []
    cos_mean = []
    cos_std = []

    for k, run in model_runs.items():
        t0 = datetime.now()
        print 'processing {}'.format(k)
        this_cos = sp.parse_STEM_var(run.aqout_path,
                                     varname='CO2_TRACER1',
                                     t0=datetime(2008, 7, 8),
                                     t1=datetime(2008, 8, 31, 23, 59, 59))
        t_data = pd.DatetimeIndex(this_cos['t'], freq='1H')

        this_t, this_mean = daily_window_stats(t_data,
                                               this_cos['data'],
                                               is_midday,
                                               np.mean)
        this_t, this_std = daily_window_stats(t_data,
                                              this_cos['data'],
                                              is_midday,
                                              np.std)
        t.append(this_t)
        cos_mean.append(this_mean)
        cos_std.append(this_std)
        print 'finished {}, time: {}'.format(k, str(datetime.now() - t0))

    t_dict = {key: value for key, value in zip(model_runs.keys(), t)}
    cos_mean_dict = {key: value for key, value in zip(model_runs.keys(),
                                                      cos_mean)}
    cos_std_dict = {key: value for key, value in zip(model_runs.keys(),
                                                     cos_std)}
    all_data_dict = {'t': t_dict,
                     'cos_mean': cos_mean_dict,
                     'cos_std': cos_std_dict}
    outfile = open(pickle_fname, 'wb')
    cPickle.dump(all_data_dict, outfile, protocol=2)
    outfile.close()

    return(all_data_dict)


def load_aqout_data(fname='/home/thilton/Data/STEM/aq_out_data.cpickle'):
    """
    load a cpickle data file containing aqout data written by assemble_data
    """
    warnings.warn('load_aqout_data is deprecated and will be removed in'
                  ' the future.  Please use '
                  'aqout_postprocess.aqout_container instead')

    import cPickle
    f = open(fname, 'rb')
    all_data = cPickle.load(f)
    f.close()
    return(all_data)


def combine_aqout_netcdf_files(fnames, fname_out="out.nc"):
    """combine multiple aqout_container netcdf files to a single
    netcdf file.

    Combines multiple aqout_container netcdf files produced by
    aqout_container.stats_to_netcdf() into a single netcdf file.  The
    datasets from each separate netcdf file are placed into netcdf
    groups.

    Depends on ncecat from the `NCO operators <http://nco.sourceforge.net/>`_

    ARGS
    fnames (list): list of strings containing full paths of netcdf
        files to be combined.
    fname_out (string): full path to the combined netcdf file to
        create.  Default is "./out.nc"

    RETURNS:
    0 on success
    """

    for f in fnames:
        if not os.path.exists(f):
            raise IOError("{} not found.".format(f))

    cmd = "ncecat --gag {} {}".format(" ".join(fnames), fname_out)
    subprocess.call(cmd, shell=True)
