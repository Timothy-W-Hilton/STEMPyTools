"""
short script to create pseudodata input.dat files for the enhanced
midwestern USA OCS fluxes and the enhanced western USA OCS fluxes
using inputdattools.  This creates the input.dat files from the
forward runs' AQOUT files, rather than their t_obs_pred.dat files.
This approach therefore allows for pseudodata at multiple timestamps,
rather than only the last timestamp that t_obs_pred.dat makes
available.
"""

import os
import os.path
import numpy as np

import inputdattools as idt

#-----
# set up STEM pseudo-observaitons
# define the dimensions of the STEM grid
nt = 504
nx = 124
ny = 124
nz = 22
# put pseudo-observations at STEM grid z = 7 (roughly 1000 m)
obs_z = 7
#take ten pseudo-observations, evenly spaced from hour 10 to hour 500
#of the three-week (504-hour) forward run.
obs_t = np.int_(np.linspace(10, 500, num=10).round())
# create a "mask" for the forward run output.  indices containing True
# will not become psueo-observations (they will be masked). Indices
# containing False will become pseudo-observations.
m = np.ones((nt,nz,nx,ny), dtype=bool);
m[obs_t, obs_z, ...] = False

wUSA_aq = os.path.join(os.getenv('STEM'),
                       'run.TWH_fwd_large_slab',
                       'output_wUSAx100',
                       'AQOUT-124x124-22levs-casa-cos_2008_2009.nc')
mwUSA_aq = wUSA_aq = os.path.join(os.getenv('STEM'),
                       'run.TWH_fwd_large_slab',
                       'output_mwUSAx1.5',
                       'AQOUT-124x124-22levs-casa-cos_2008_2009.nc')

print 'parsing western USA x 100 AQOUT'
wUSA_id = idt.STEMInputDat.cos_from_aqout(wUSA_aq, m)
print 'writing western USA x 100 AQOUT input.dat'
wUSA_id.write_file(fname='wUSAx100_10timesteps.input.dat',
                   fdir=os.path.join(os.getenv('HOME'),
                                     'Data', 'STEM', 'input'))

print 'parsing midwestern USA x 1.5 AQOUT'
mwUSA_id = idt.STEMInputDat.cos_from_aqout(mwUSA_aq, m)
print 'writing midwestern USA x 1.5 input.dat'
mwUSA_id.write_file(fname='mwUSAx1.5_10timesteps.input.dat',
                   fdir=os.path.join(os.getenv('HOME'),
                                     'Data', 'STEM', 'input'))
