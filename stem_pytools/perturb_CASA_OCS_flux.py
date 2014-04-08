import os
import os.path
import sys
import shutil
import netCDF4
import numpy as np
import numpy.ma as ma

import STEM_vis
import STEM_parsers

class customSTEMFlux(object):
    """
    Class to parse in, manipulate, and write out a STEM OCS flux
    netcdf file.  This is intended to facilitate dialing the STEM flux
    up or down for sub regions of North America in order to support
    twin experiment studies.
    """

    def __init__(self,
                 fname_flux_orig,
                 fname_out,
                 manip_mask,
                 manip_factor):
        self.fname_flux_orig = fname_flux_orig
        self.fname_out = fname_out
        self.manip_mask = manip_mask
        self.manip_factor = manip_factor

    def manipulate(self, copy=True):
        """
        Apply the requested manipulation to the specified file:
        multiply the fluxes in fname_flux_orig at the locations where
        manip_mask is true by manip_factor
        """
        if copy:
            shutil.copyfile(self.fname_flux_orig, self.fname_out)
        nc = netCDF4.Dataset(self.fname_out, 'a', format='NETCDF4')
        flx = nc.variables['cos']

        new_flx = flx[:]
        new_flx[..., self.manip_mask] *= self.manip_factor
        nc.variables['cos'][:] = new_flx

        nc.close()
        return(new_flx)

def make_midwest_mask(fname_coords):
    """
    Define a mask for the 124x124 STEM grid containing True for a
    rectangle covering (roughly) Iowa and the northern half of
    Illinois, and containing False everywhere else.
    """
    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(fname_coords)
    lat_range = np.array((40,44))  #deg N
    lon_range = np.array((-96, -87)) #deg W

    lat_mask = ma.masked_inside(lat, *lat_range)
    lon_mask = ma.masked_inside(lon, *lon_range)
    mask = np.logical_and(lat_mask.mask, lon_mask.mask)

    return(mask)

def make_western_mask(fname_coords, fname_top):
    """
    All my STEM optimizations want the emissions scaling factor to be
    uniformly 1.0 almost everywhere west of (roughly) the front range
    of the Rockies.  Given a t_obs_pred.dat file from such a STEM
    optimization, define a mask which is true in these locations (both
    West of Rockies and where emi_fac is 1.0) and false elsewhere.
    """
    top = STEM_parsers.parse_tobspred(fname_top)
    top = STEM_vis.grid_tobspred_data(top, which_data='emi_fac')

    lon, lat, topo = STEM_parsers.parse_STEM_coordinates(fname_coords)

    mask = np.logical_and(abs(top - 1.0) < 0.001, lon < -95)
    # mask out this little blip in prairie that slips through
    mask[np.logical_and(lat > 40, lon > -100)] = False
    #mask = lon < -110

    return(mask)

if __name__ == '__main__':
    mask_mw = ptb.make_midwest_mask(os.path.join(os.getenv('SARIKA_INPUT'),
                                                 'TOPO-124x124.nc'))
    mask_W = ptb.make_western_mask(
        os.path.join(os.getenv('SARIKA_INPUT'),
                     'TOPO-124x124.nc'),
                     os.path.join('/home', 'thilton',
                                  'Stem_emi2_onespecies_big_ocssib',
                                  'run.TWH_opt_LargeSlab_StrongPrior_1.5xFlux_21day',
                                  't_obs_pred_013.dat'))
    
    original_CASA_fluxes = os.path.join(os.getenv('SARIKA_INPUT'),
                                        'surfem-124x124-casa-cos_2008_2009.nc')
    
    mw = ptb.customSTEMFlux(
        original_CASA_fluxes,
        os.path.join(os.getenv('HOME'),
                     'Data', 'STEM', 'input',
                     'surfem-124x124-casa-cos_2008_2009_mwUSAx1.5.nc'),
                     mask_mw,
                     1.5)
    flx_mw = mw.manipulate()
    
    w = ptb.customSTEMFlux(
        original_CASA_fluxes,
        os.path.join(os.getenv('HOME'),
                     'Data', 'STEM', 'input',
                     'surfem-124x124-casa-cos_2008_2009_wUSAx100.nc'),
                     mask_W,
                     100.0)
    flx_w = w.manipulate()
