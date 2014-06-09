"""A collection of functions to parse STEM input and output files"""

from __future__ import division
import os.path
from glob import glob
import re
import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime
from netCDF4 import Dataset

#--------------------------------------------------
# helper classes, functions

class BlockReader(object):
    """
    class to read ASCII text in n-line "blocks"
    """
    def __init__(self, f, n=1):
        self.f = f
        self.n = n
    def __iter__(self):
        return self
    def next(self):
        return [self.f.next() for i in xrange(self.n)]

def string_strip(lst):
    """
    apply strip to each string in a list of strings
    """
    return([x.strip() for x in lst])


def file_len(fname):
    """returns the number of lines contained in the file fname."""
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class StemRunFile(object):
    """
    Class to parse variable/value pairs from STEM run files.
    CLASS ATTRIBUTES:
    fname: the full path to the STEM run file
    lines: list of str; all lines from the run file containing
        variable assignments
    vars: dict; variable name, value pairs
    t_start: datetime.datetime; the starting time of the STEM simulation
    t_end: datetime.datetime; the ending time of the STEM simulation
    """
    def __init__(self, fname):
        """
        Class constructor.  Parses the specified file and populates
        lines and vars.

        INPUTS
        fname: str; full path to a STEM run file
        """
        self.fname=fname
        self.parse_to_list()
        self.trim_lines()
        self.create_dict()
        self.sub_vars()
        self.calc_run_start()
        self.calc_run_end()
    def parse_to_list(self):
        f = open(self.fname)
        self.lines = f.readlines()
        f.close()
    def trim_lines(self):
        """
        discard lines from the run file that do not contain
        variable/value assigments.  Lines that do not begin with a '#'
        and do contain '=' or 'setenv' are assumed to be assignments
        """
        #get rid of comments
        self.lines = [ln for ln in self.lines if ln[0] is not '#']
        #keep variable assigments (shell and environment)
        self.lines = [ln for ln in self.lines if
                      ('=' in ln) or ('setenv' in ln)]
    def create_dict(self):
        """
        for each line in the format
        "var=val"
        or
        "sentenv var val"
        put var and val into a dict
        """
        #re to match "var = val", with arbitrary whitespace around the =
        m_eq = [re.search('(?P<var>\w+)=(?P<val>.+)', ln) for
                ln in self.lines]
        m_eq = [mat for mat in m_eq if mat is not None]
        #re to match "setenv var val", with arbitrary whitespace separating
        m_env = [re.search('setenv\s(?P<var>.+)\s(?P<val>.+)', ln) for
                 ln in self.lines]
        m_env = [mat for mat in m_env if mat is not None]

        #combine the two lists of dicts into one big dict
        merged_dict = {}
        for m in m_eq:
            d = m.groupdict()
            merged_dict[d['var']] = d['val']
        for m in m_env:
            d = m.groupdict()
            merged_dict[d['var']] = d['val']
        self.vars=merged_dict
    def sub_vars(self):
        """
        substitute environment variables referenced in the run file
        with their values when the specify paths.  For variables in
        the format $ABC_DEF, try first to replace from the other
        variables defined in the run file.  If the variable is not
        present, look within os environment variables.  If both fail
        leave the variable unchanged.
        """
        for k in self.vars.keys():
            # find environment variables
            match = re.search('\$(?P<varname>[A-Z0-9_]+)', self.vars[k])
            if match is not None:
                varname = match.group('varname')
                if ((varname in self.vars.keys()) and
                    (os.path.exists(self.vars[varname]))):
                    full_val = self.vars[varname]
                    self.vars[k] = self.vars[k].replace('$' + varname,
                                                        full_val)
                elif ((os.getenv(varname) is not None) and
                    (os.path.exists(os.getenv(varname)))):
                    full_val = os.getenv(varname)
                    self.vars[k] = self.vars[k].replace('$' + varname,
                                                        full_val)
    def calc_run_start(self):
        t = (datetime.datetime.strptime(self.vars['istday'], '%Y %m %d') +
              datetime.timedelta(hours=int(self.vars['isthr'])))
        self.t_start = t
    def calc_run_end(self):
        dt = datetime.timedelta(hours=int(self.vars['iperiod']))
        self.t_end = self.t_start + dt


#--------------------------------------------------
# parser functions

def parse_inputdat(fname):
    """parse the STEM input.dat file specified by fname to a pandas
    data frame.

    RETURNS: a data frame with the variables t, x, y, z, and COS.
    """
    input_dat = pd.read_csv(fname,
                            sep=None,
                            header=None,
                            skiprows=4,
                            names=('t', 'x', 'y', 'z', 'COS'))
    return(input_dat)

def parse_tobspred(fname, inputdat_fname=None, emi_fac_l_no=None):
    """ parse model and observed OCS concentrations and emissions
    scaling factors ('emi_fac') from a STEM t_obs_pred.dat file
    specified by fname to a two pandas data frames.
    PARAMETERS
    ----------
    fname; string: full path to the t_obs_pred.dat file to be read
    inputdat_fname; string, optional: full path to the input.dat file
       that drove the generation o f the t_obs_pred.dat file.  If
       specified the input.dat file is parsed to determine how many
       lines of concentrations there are in the t_obs_pred.dat file.
    emi_fac_l_no: integer, optional: The line number of the first line
       in t_obs_pred.dat that contains emissions scaling factors.
       Ignored if inputdat_fname is specified.  Must be specified if
       inputdat_fname is not specified.

    RETURNS: a dict with two entries: ocs_conc and emi_fac.  Each
    entry contains a pandas data frame.  ocs_conc contains time stamps
    (t), observed OCS concentrations (obs), and model OCS
    concentrations (mod).  emi_fac contains STEM grid X coordinate
    (x), STEM grid Y coordinate (y), and the emissions scaling factor
    (emi_fac).

    Note: The input.dat file is used to identify the line in the
    t_obs_pred.dat file where the data change from concentrations to
    emissions scaling factors.  The only way I could think of to do
    this completely from the information within the t_obs_pred.dat
    file is to identify the first line where the first number is an
    integer, rather than a floating point number in exponential
    notation.  This approach would be vulnerable to a change in output
    format within STEM, though; therefore I decided to go with using
    input.dat.
    """
    if inputdat_fname is not None:
        n_hdr_lines = 4;
        emi_fac_l_no = file_len(inputdat_fname) - n_hdr_lines
    elif emi_fac_l_no is None:
        raise TypeError('Neither inputdat_fname nor'
                        'emi_fac_l_no were specified')

    ocs_conc = pd.read_csv(fname,
                           sep='[\s]*',
                           header=None,
                           skipinitialspace=False,
                           nrows=emi_fac_l_no,
                           names=('t', 'obs', 'mod'))

    emi_fac = pd.read_csv(fname,
                          sep='[\s]*',
                          header=None,
                          skipinitialspace=False,
                          skiprows=emi_fac_l_no,
                          names=('x', 'y', 'emi_fac'))
    return({"ocs_conc" : ocs_conc,
            "emi_fac" : emi_fac})

def parse_reportopt(fname, block_sz=11):
    """
    parse a STEM inverse run Report.opt output file to a pandas data frame.

    RETURNS: a pandas data frame with variables:
       it: model iteration
       mod_runs: number of model runs
       cost: the value of the cost functions
       misfit: the 'misfit' component of the cost function
       bckg: the 'background' component of the cost function
       task: the L-BFGS 'task' for the iteration
       """
    block_list = []

    with open(fname, 'r') as f:
        for block in BlockReader(f, 11):
            # build a regex to match floating point numbers
            re_flt = r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
            for this_line in block:
                if re.search("LBFGS-TASK", this_line):
                    task=re.split('LBFGS-TASK =', this_line)[1].strip()
                if re.search("It=", this_line):
                    it, nruns = string_strip(re.findall(re_flt, this_line))
                if re.search("Cost", this_line):
                    cost, = string_strip(re.findall(re_flt, this_line))
                if re.search("Misfit", this_line):
                    misfit, bckg = string_strip(re.findall(re_flt, this_line))
            block_list.append({'task':task,
                               'it':int(it),
                               'nruns':int(nruns),
                               'cost':float(cost),
                               'misfit':float(misfit),
                               'bckg':float(bckg)})
    f.close()
    df = pd.DataFrame(block_list)
    return(df)

def get_all_tobspred_fnames(run_dir):
    """returns list of full paths to all files in the specified
    directory matching t_obs_pred*.dat.  The results are sorted
    lexically"""
    file_list = glob(os.path.join(run_dir,
                                  't_obs_pred*.dat'))
    file_list = sorted(file_list)
    return(file_list)

def parse_all_emifac(run_dir, mask_ones=True):
    """parse all emi_fac values from all t_obs_pred_ABC.dat files
    present within a specified directory into a numpy array, one
    column per STEM iteration.  Provides an option (on by default) to
    mask values emi_fac values equal to 1.0."""
    emifac = None
    fnames = get_all_tobspred_fnames(run_dir)
    if fnames:
        emifac_list = [parse_tobspred(f) for f in fnames]
        emifac = [x['emi_fac']['emi_fac'].values for x in emifac_list]
        emifac = np.transpose(np.array(emifac))
        #mask emifac values == 1.0
        if mask_ones:
            emifac = np.ma.masked_array(emifac, (np.abs(emifac) - 1.0) < 1e-10)
    return(emifac)

def parse_STEM_tflag(nc_fname, out_format='datetime'):
    """
    parse the TFLAG variable of a STEM IO/API file to
    datetime.datetime values.  TFLAG is the time variable in STEM
    input and output files.  It's format is a text string:
    YYYYDDD,HHMMSS.  The specified file must contain the variable
    TFLAG.

    PARAMETERS
    ----------
    nc_fname: string; the full path to the IO/API file.
    out_format: {datetime}|hour: format to return time.
    """
    SECONDS_PER_HOUR = 60*60
    try:
        nc = Dataset(nc_fname, 'r', format='NETCDF4')
    except:
        print('error opening {}'.format(nc_fname))
        raise
    # read timestamps to datetime.datetime
    t = np.squeeze(nc.variables['TFLAG'])
    t_dt = np.array(
        ([datetime.datetime.strptime(str(this[0]) +
                                     str(this[1]).zfill(6), '%Y%j%H%M%S')
                                     for this in t]))
    nc.close()

    if out_format is 'hour':
        # conver the datetime objects to hours past the first timestamp.
        t0 = t_dt[0]
        # t_dt has dtype 'O'; numpy refuses to convert these to floats
        # after hour calculation, so initialize a new array of dtype
        # float to hold hours.
        t_hr = np.empty_like(t_dt, dtype=float)
        for i in range(t_dt.size):
            td = t_dt[i] - t0
            t_hr[i] = td.total_seconds() / SECONDS_PER_HOUR
        t_dt = t_hr
    return(t_dt)

def parse_STEM_var(nc_fname=None,
                   t_idx=None,
                   z_idx=None,
                   t0=None,
                   t1=None,
                   varname=None):
    """ Parse a STEM variable from a STEM I/O API netcdf file.
    varname (type str) must be a variable in the netcdf file. The file
    must also contain a variable TFLAG containing timestamps in the
    format <YYYYDDD,HHMMSS>.

    There are two ways so specify the time slices to extract:
    (1) specify t_idx: array-like indices directly into the time
    dimension of the IO/API file.  If t_idx is specified t1 and t0 are
    ignored.
    (2) specify, t0 and/or t1: (datetime.datetime) restrict the
    returned data to timestamps that satisfy t0 <= timestamp <= t1.
    If none of t_idx, t0,and t1 are not specified then all timestamps
    are returned.  If t_idx is specified t1 and t0 are ignored.

    z_idx specifies the vertical layers to extract.  Must be None (all
    vertical layers) , a scalaer integer, or a tuple.

    RETURNS a dict with keys 'data' and 't'.  'data' contains the
    values in varname (np.ndarray) and 't' contains the timestamps
    (datenum.datenum objects)."""
    try:
        nc = Dataset(nc_fname, 'r', format='NETCDF4')
    except:
        print('error opening {}'.format(nc_fname))
        raise
    t_dt = parse_STEM_tflag(nc_fname)
    if t_idx is None:
        # find the requested timestamps
        if t0 is None:
            t0 = t_dt.min()
        if t1 is None:
            t1 = t_dt.max()
        t_idx = (t_dt >= t0) & (t_dt <= t1)
    if z_idx is None:
        z_idx = np.arange(nc.variables[varname].shape[1])
    elif type(z_idx) is not tuple:
        z_idx = (z_idx,)
    # retrieve the requested [OCS] data
    data = nc.variables[varname][t_idx, z_idx, :, : ]
    nc.close()
    return({'data':data, 't':t_dt[t_idx]})

def parse_STEM_coordinates(topo_fname):
    """Parse STEM grid latitude and longitude."""
    try:
        topo = Dataset(topo_fname, 'r', format='NETCDF4')
    except:
        print('error opening {}'.format(topo_fname))
        raise

    lat = np.squeeze(topo.variables['LAT'])
    lon = np.squeeze(topo.variables['LON'])
    topo = np.squeeze(topo.variables['TOPO'])

    return(lon, lat, topo)

def get_CO2grossflux_varname(nc_fname):
    """
    determine whether the CO2 gross flux variable is 'GPP' or 'GEE'.

    examine variables in nc_fname and return 'GEE' if present.  If GEE
    is not present, return 'GPP' if present.  If neither 'GEE' or
    'GPP' are present return None.
    """
    try:
        nc = Dataset(nc_fname)
    except:
        print('error opening {}'.format(nc_fname))
        raise

    if 'GEE' in nc.variables.keys():
        return('GEE')
    elif 'GPP' in nc.variables.keys():
        return('GPP')
    else:
        return(None)
