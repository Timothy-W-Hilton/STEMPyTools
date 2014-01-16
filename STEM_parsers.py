### functions to parse STEM input and output files

from __future__ import division
import os.path
import pandas as pd

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def parse_inputdat(fname):
    """parse specified STEM input.dat file to a pandas data frame"""
    input_dat = pd.read_csv(fname,
                            sep=None,
                            header=None,
                            skiprows=4,
                            names=('t', 'x', 'y', 'z', 'COS'))
    return(input_dat)

def parse_tobspred(fname):
    """ parse model and observed OCS concentrations and emissions
    scaling factors ('emi_fac') from a specified t_obs_pred.dat file"""

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
