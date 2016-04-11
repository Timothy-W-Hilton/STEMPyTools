import numpy as np
import numpy.ma as ma
import os
import os.path

import domain


def calc_STEM_COS_drawdown(aqout_conc,
                           topo_fname=None,
                           wrfheight_fname=None,
                           lo_height_agl=2000,
                           hi_height_agl=4000):
    """
    calculate STEM COS vertical drawdown from AQOUT.  As per
    conversation with Elliott on 30 Oct 2014, this function defines
    vertical drawdown as (mean [COS], 0 m to 2000 m) minus (mean [COS]
    > 4000 m).

    INPUTS
    aqout_conc: ndarray; AQOUT [COS], molecules m-3
    topo_fname: full path to STEM topo file.  Default is
       $SARIKA_INPUT/TOPO-124x124.nc.
    wrfheight_fname: full path to STEM WRF height file.  Default is
       $SARIKA_INPUT/wrfheight-124x124-22levs.nc
    lo_height_agl: height above ground level, in meters, of the top of
       the "surface" altitdude bin.  Default is 2000 m.
    hi_height_agl: height above ground level, in meters, of the bottom of
       the "high" altitdude bin.  Default is 4000 m.
    """
    MCLS_M3_2_PPTV = 1e12

    if topo_fname is None:
        topo_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                  'TOPO-124x124.nc')
    if wrfheight_fname is None:
        wrfheight_fname = os.path.join(os.getenv('SARIKA_INPUT'),
                                       'wrfheight-124x124-22levs.nc')
    d = domain.STEM_Domain(fname_topo=topo_fname)
    d.get_STEMZ_height(wrfheight_fname)

    # tile agl a 3D array. Tile it out to four dimensions so it has
    # the same number of time stamps as aqout_conc.
    agl = np.tile(d.agl, (aqout_conc.shape[0], 1, 1, 1))

    stem_z = 1  # second axis of AQOUT is z level

    lo_cos = ma.masked_where(agl >= lo_height_agl,
                             aqout_conc).mean(axis=stem_z)
    hi_cos = ma.masked_where(agl <= hi_height_agl,
                             aqout_conc).mean(axis=stem_z)

    dd = (hi_cos - lo_cos) * MCLS_M3_2_PPTV
    dd = dd[:, np.newaxis, ...]

    abl_mean = lo_cos * MCLS_M3_2_PPTV
    free_trop_cos = hi_cos * MCLS_M3_2_PPTV

    return(dd)
