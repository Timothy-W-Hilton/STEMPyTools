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

def parse_tobspred(fname):
    """ parse model and observed OCS concentrations and emissions
    scaling factors ('emi_fac') from a STEM t_obs_pred.dat file
    specified by fname to a two pandas data frames.

    RETURNS: a dict with two entries: ocs_conc and emi_fac.  Each
    entry contains a pandas data frame.  ocs_conc contains time stamps
    (t), observed OCS concentrations (obs), and model OCS
    concentrations (mod).  emi_fac contains STEM grid X coordinate
    (x), STEM grid Y coordinate (y), and the emissions scaling factor
    (emi_fac).
    """

    nlines = file_len(fname)
    if (nlines % 2) != 0:
        raise "t_obs_pred.dat must have even number of lines!"

    ocs_conc = pd.read_csv(fname,
                           sep='[\s]*',
                           header=None,
                           skipinitialspace=False,
                           nrows=int(nlines / 2),
                           names=('t', 'obs', 'mod'))

    emi_fac = pd.read_csv(fname,
                          sep='[\s]*',
                          header=None,
                          skipinitialspace=False,
                          skiprows=int(nlines / 2),
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

def parse_STEM_var(nc_fname=None, t0=None, t1=None, varname=None):
    """ Parse a STEM variable from a STEM I/O API netcdf file.
    varname (type str) must be a variable in the netcdf file. The file
    must also contain a variable TFLAG containing timestamps in the
    format <YYYYDDD,HHMMSS>.  If specified, t0 and t1
    (datetime.datetime) restrict the returned data to timestamps that
    satisfy t0 <= timestamp <= t1.  If t0 and t1 are not specified
    then all timestamps are returned.

    RETURNS a dict with keys 'data' and 't'.  'data' contains the
    values in varname (np.ndarray) and 't' contains the timestamps
    (datenum.datenum objects)."""
    nc = Dataset(nc_fname, 'r', format='NETCDF4')
    # read timestamps to datetime.datetime
    t = np.squeeze(nc.variables['TFLAG'])
    t_dt = np.array(
        ([datetime.datetime.strptime(str(this[0]) +
                                     str(this[1]).zfill(6), '%Y%j%H%M%S')
                                     for this in t]))
    # find the requested timestamps
    t_idx = (t_dt >= t0) & (t_dt <= t1)
    # retrieve the requested [OCS] data
    data = nc.variables[varname][t_idx, 0, :, : ]
    return({'data':data, 't':t_dt[t_idx]})

def parse_STEM_coordinates(topo_fname):
    """Parse STEM grid latitude and longitude."""
    topo = Dataset(topo_fname, 'r', format='NETCDF4')

    lat = np.squeeze(topo.variables['LAT'])
    lon = np.squeeze(topo.variables['LON'])
    topo = np.squeeze(topo.variables['TOPO'])

    return(lon, lat, topo)
