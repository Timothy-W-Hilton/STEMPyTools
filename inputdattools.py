"""
inputdattools: STEM input.dat tools
A python module to provide tools to work with STEM input.dat files.
Functionality includes:
(1) assemble pseudo-observation input.dat from a STEM forward run
    AQOUT IO/API file
(2) assemble input.dat from OCS observations

Provides classes:
(1) STEMInputDat
"""

import os, os.path
import numpy as np
import numpy.ma as ma
import netCDF4
import pandas as pd

import STEM_parsers

class STEMInputDat(object):
    """
    Container class to hold the information needed to write a STEM
    input.dat file: species concentration observations and
    accompanying hour and x, y, and z indices.  Currently only
    supports one species.
    """
    def __init__(self,
                 hr=None,
                 t=None,
                 x=None,
                 y=None,
                 z=None,
                 obs=None,
                 n_total_spc=1,
                 hdr_note='##t, x, y, z, COS (ppbv)',
                 uncrt=(0.8,),
                 spc_idx=1):
        """
        Class constructor.
        PARAMETERS
        ----------
        x, y, z: N-element np.ndarray:  x, y, z indices of the observations
        t:  N-element np.ndarray; STEM timestep of the observations
        hr: N-element np.ndarray; hour of model run corresponding to
            observations
        obs: N by n_total_spc np.ndarray; columns species concentrations
        n_total_spc: int; number of species in the input.dat file.
        hdr_note: comment to appear in the input.dat file header
        uncrt: tuple of ints, length n_total_spc: the uncertainty corresponding
            to each species' observations.
        spc_idx: species index that matches aq_property.dat
        """
        self.n_total_spc = n_total_spc
        self.hdr_note = hdr_note
        self.uncrt = uncrt
        self.spc_idx = spc_idx
        self.t = t
        self.hr = hr
        self.x = x
        self.y = y
        self.z = z
        self.obs = obs

        # TWH design notes
        # if masked, (1) use np.nonzero to get indices where mask is false
        #            (2) make DataFrame of indices, value
        #            (3) write DataFrame to input.dat.  Probably want method
        #            to write input.dat based on header info and a data frame
        # if ndarray: (1) make it a masked array where the mask is everywhere false
        #             (2) proceed as above

    @classmethod
    def cos_from_inputdat(cls, inputdat_fname):
        """
        parse an input.dat file to a STEMInputDat object
        """

        f = open(inputdat_fname, 'r')
        n_obs, n_total_spc = [int(i) for i in f.readline().strip().split()]
        hdr_note = f.readline().strip()
        uncrt, spc_idx = f.readline().strip().split()
        uncrt = float(uncrt)
        spc_idx = int(spc_idx)
        f.close()

        input_dat = pd.read_csv(inputdat_fname,
                                sep=None,
                                header=None,
                                skiprows=4,
                                names=('t', 'x', 'y', 'z', 'COS'))

        obj = cls(t=input_dat['t'].get_values(),
                  x=input_dat['x'].get_values(),
                  y=input_dat['y'].get_values(),
                  z=input_dat['z'].get_values(),
                  obs=input_dat['COS'].get_values(),
                  n_total_spc=n_total_spc,
                  hdr_note=hdr_note,
                  uncrt=uncrt,
                  spc_idx=spc_idx)
        return(obj)

    @classmethod
    def cos_from_aqout(cls,
                       aq_fname,
                       aq_mask=None):
        """
        Method to build a STEM input.dat containing pseudo-observations
        from a subset of a STEM forward run AQOUT IO/API file.
        PARAMETERS
        ----------
        aq_fname: full path to the AQOUT file containing the pseudo-observations.
        aq_mask: ndarray, dtype bool; array of the same shape as the
           observations in AQOUT.  The OCS concentrations in AQOUT
           where aq_mask is True will become pseudo-observations.  The
           rest will be ignored.  If aq_mask is not specified all
           values from AQOUT will be kept.
        """
        # factor to convert AQOUT [OCS] units of molecules OCS cm^-3 to ppbv
        unit_convert_factor = 1e9
        t_hr = np.int_(STEM_parsers.parse_STEM_tflag(aq_fname,
                                                     out_format='hour'))

        nc = netCDF4.Dataset(aq_fname, 'r')
        nc_ocs = nc.variables['CO2_TRACER1']

        if aq_mask is not None:
            # find observations that are not masked
            t_idx, z_idx, x_idx, y_idx = np.nonzero(np.logical_not(aq_mask))
        else:
            t_idx, z_idx, x_idx, y_idx = map(np.range, nc_ocs.shape)
        vals = np.empty_like(t_idx, dtype=float)
        vals[:] = np.NaN

        for i in range(t_idx.size):
            vals[i] = (nc_ocs[t_idx[i], z_idx[i], x_idx[i], y_idx[i]] *
                       unit_convert_factor)


        nc.close()

        obj = cls(obs=vals,
                  hr=t_hr[t_idx],
                  x=x_idx + 1,  # input.dat uses one-based indexing
                  y=y_idx + 1,  # input.dat uses one-based indexing
                  z=z_idx + 1,  # input.dat uses one-based indexing
                  hdr_note='##t, x, y, z, COS (ppbv) | pseudo-observations')
        return(obj)

    def write_file(self,
                   fname='input.dat',
                   fdir=None):
        """
        write an input.dat file to disk from the contents of a
        STEMInputDat object.
        parameters:
        """
        hdr_str = ("{n_obs}  {n_total_spc}\n"
                   "{hdr_note}\n"
                   "{spc_uncrt} {n_total_spc}\n"
                   "{n_spc}".format(n_obs=self.hr.size,
                                      n_total_spc=self.n_total_spc,
                                      hdr_note=self.hdr_note,
                                      spc_uncrt=self.uncrt[0],
                                      n_spc=self.spc_idx))
        arr = np.concatenate((self.hr[...,np.newaxis],
                              self.x[...,np.newaxis],
                              self.y[...,np.newaxis],
                              self.z[...,np.newaxis],
                              self.obs[...,np.newaxis]), axis=1)
        if fdir is not None:
            fname = os.path.join(fdir, fname)
        np.savetxt(fname,
                   arr,
                   fmt=('%d', '%d', '%d', '%d', '%0.10f'),
                   delimiter='\t',
                   header=hdr_str,
                   comments='')


# @classmethod
# def cos_from_aqout(cls,
#                    aq_fname,
#                    aq_mask=None):
#     """
#     Method to build a STEM input.dat containing pseudo-observations
#     from a subset of a STEM forward run AQOUT IO/API file.
#     INPUTS
#     aq_fname: full path to the AQOUT file containing the pseudo-observations.
#     aq_mask: ndarray, dtype bool; array of the same shape as the
#     observations in AQOUT.  The OCS concentrations in AQOUT where
#     aq_mask is True will become pseudo-observations.  The rest will
#     be ignored.  If aq_mask is not specified all values from AQOUT
#     will be kept.
#     """
#     if aq_mask is not None:
#         # find observations that are not masked
#         t_idx, z_idx, x_idx, y_idx = np.nonzero(np.logical_not(aq_mask))
#     else:
#         #keep all t, Z values
#         t_idx = None
#         z_idx = None

#     import pdb; pdb.set_trace()
#     # parse the observations from the IO/API file
#     # netCDF4 currentlydoes not allow multi-dimensional indexing:
#     #     https://code.google.com/p/netcdf4-python/issues/detail?id=106
#     aq = STEM_parsers.parse_STEM_var(aq_fname,
#                                      t_idx=tuple(t_idx),
#                                      z_idx=tuple(z_idx),
#                                      varname='CO2_TRACER1')
#     import pdb; pdb.set_trace()
#     if aq_mask is None:
#         t_idx, z_idx, x_idx, y_idx = np.unravel_index(indices=np.arange(aq['data'].size), dims=np.arange(aq['data'].shape))
#         ocs = aq['data'][t_idx, z_idx, x_idx, y_idx]

#     #return((t_idx, z_idx, x_idx, y_idx, ocs))
#     # the DataFrame gives a Memory error (probably because it
#     # copies the indices in memory.  the indices have to be int64s
#     # to span the size of a 21-day input array, so 5 100K+ element
#     # arrays of 64 bit values is 630MB each...  times five and
#     # copied I guess it's not surprising that memory runs out
#     data = pd.DataFrame({'t':t_idx,
#                          'Z':z_idx,
#                          'X':x_idx,
#                          'Y':y_idx,
#                          'OCS':aq['data'][t_idx, z_idx, x_idx, y_idx]})

#     # obj = cls(obs=data,
#     #           hdr_note='##t, x, y, z, COS (ppbv) | pseudo-observations')
#     # return(obj)
