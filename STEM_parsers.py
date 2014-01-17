"""A collection of functions to parse STEM input and output files"""

from __future__ import division
import os.path
import re
import pandas as pd

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

