import os
import os.path
import numpy as np
import pandas as pd
import itertools
from datetime import datetime, timedelta
import cPickle

from stem_pytools import ecampbell300_data_paths as edp
from stem_pytools import STEM_parsers as sp


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
        is the output of stem_pytools.ecampbell300_data_paths.get_runs()
    pickle_fname: full path of the cpickle file to create.  If
        unspecified the default is
        /home/thilton/Data/STEM/aq_out_data.cpickle

    OUTPUT PARAMETERS:
    all_data: dict containing (1), (2), and (3) above.

    """
    if model_runs is None:
        model_runs = edp.get_runs()

    t = []
    cos_mean = []
    cos_std = []

    for k, run in model_runs.items():
        t0 = datetime.now()
        print 'processing {}'.format(k)
        this_cos = sp.parse_STEM_var(run.aqout_path,
                                     varname='CO2_TRACER1',
                                     t0=datetime(2008, 7, 1),
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
    outfile = open('/home/thilton/Data/STEM/aq_out_data_BASC.cpickle', 'wb')
    cPickle.dump(all_data_dict, outfile, protocol=2)
    outfile.close()

    return(all_data_dict)


def load_aqout_data(fname='/home/thilton/Data/STEM/aq_out_data.cpickle'):
    """
    load a cpickle data file containing aqout data written by assemble_data
    """
    import cPickle
    f = open(fname, 'rb')
    all_data = cPickle.load(f)
    f.close()
    return(all_data)


def demo():
    fname = os.path.join('/', 'home', 'thilton',
                         'Stem_emi2_onespecies_big_ocssib',
                         'STEM_Runs_LRU_Paper',
                         'CanIBIS_LRU1.61',
                         'output',
                         'AQOUT-124x124-22levs-CanIBIS_fCOS_LRU1.61.nc')

    fname = os.path.join('/', 'Users', 'tim', 'work', 'Data', 'STEM',
                         'output',
                         'AQOUT-124x124-22levs-CanIBIS_fCOS_LRU1.61.nc')

    t_ex = [datetime.now()]
    cos = sp.parse_STEM_var(fname, varname='CO2_TRACER1')
    t_ex.append(datetime.now())
    t_mn, cos_md_mean = daily_window_stats(
        pd.DatetimeIndex(cos['t'], freq='1H'),
        cos['data'],
        is_midday,
        np.mean)
    t_ex.append(datetime.now())
    t_mn, cos_md_std = daily_window_stats(
        pd.DatetimeIndex(cos['t'], freq='1H'),
        cos['data'],
        is_midday,
        np.std)
    t_ex.append(datetime.now())

    return(t_ex, t_mn, cos_md_mean, cos_md_std)
