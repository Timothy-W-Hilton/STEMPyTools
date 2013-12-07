import matplotlib as mpl
mpl.use('Agg')

import os.path
import numpy as np
import numpy.ma as ma
from datetime import datetime
from mpl_toolkits.basemap import maskoceans
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import STEM_vis
import na_map

read_data = True
if read_data:
    topo_fname = '/mnt/home10/skulkarni/StemData21Jul2013/input/TOPO-124x124.nc'
    aq_fname = os.path.join('/mnt/home10/thilton/',
                            'Stem_emi2_onespecies_big_ocssib/run.TWH/output/',
                            'AQOUT-124x124-22levs-casa-cos_2008_2009.nc')
    em_fname = os.path.join('/mnt/home10/skulkarni/StemData21Jul2013/input',
                            'surfem-124x124-casa-cos_2008_2009.nc')

    t_0 = datetime(2008, 6, 1)
    t_end = datetime(2008, 6, 30)


    [ ocs, t ] = STEM_vis.parse_STEM_var(nc_fname=aq_fname,
                                         t0=t_0,
                                         t1=t_end,
                                         varname='CO2_TRACER1')

    [ ocs_flx, t_ocs_flx ] = STEM_vis.parse_STEM_var(nc_fname=em_fname,
                                                     t0=t_0,
                                                     t1=t_end,
                                                     varname='cos')

    lon, lat = STEM_vis.parse_STEM_coordinates( topo_fname )

    mask = maskoceans( lon, lat, np.zeros( lon.shape ), inlands = False )
    mask = mask.mask  # get rid of the dummy zeros
    mask = np.tile( mask, [ ocs.shape[0], 1, 1 ] )
    ocs = ma.masked_array( ocs, mask = mask )
    ocs_flx = ma.masked_array( ocs_flx, mask = mask )

t_lab = None #initialize time label
#for i in range( ocs.shape[0] ):
for i in range( 39,42 ):
    t0 = datetime.now()
    # m = STEM_vis.initialize_STEM_map()
    # cs = m.add_ocs_contour_plot(lon,
    #                             lat,
    #                             ocs[ i, :, : ],
    #                             vmin=ocs.min(),
    #                             vmax=ocs.max(),
    #                             t_str=t[i],
    #                             cbar_t_str='mcls cm$^{-3}$')

    # fname = datetime.strftime( t[i], '/home/thilton/Plots/STEM_Maps/%Y-%m-%dT%H%M_STEM_OCS.png' )
    # print 'saving ' + fname + '(' + str(datetime.now() - t0) + ')'
    # plt.savefig( fname )
    # plt.close( m.fig )

## TODO
##    - documentation

    # m_flx = STEM_vis.initialize_STEM_map()
    # cs = m_flx.add_ocs_contour_plot(lon,
    #                                 lat,
    #                                 ocs_flx[ i, :, : ],
    #                                 vmin=ocs_flx.min(),
    #                                 vmax=ocs_flx.max(),
    #                                 t_str=t_ocs_flx[i],
    #                                 cmap=cm.get_cmap('Blues_r'),
    #                                 cbar_t_str='mol m$^{-2}$ s$^{-1}$')

    fig_combined = plt.figure(figsize=(16, 8))
    ax_conc = plt.subplot(1,2,1)
    ax_flx = plt.subplot(1,2,2)

    m_conc = na_map.NAMapFigure(map_axis=ax_conc)
    cs = m_conc.add_ocs_contour_plot(lon,
                                lat,
                                ocs[ i, :, : ],
                                vmin=ocs.min(),
                                vmax=ocs.max(),
                                t_str=t[i],
                                cbar_t_str='mcls cm$^{-3}$')
    plt.colorbar(mappable=cs,
                 ax=ax_conc,
                 format='%0.2e')

    m_flx = na_map.NAMapFigure(map_axis=ax_flx)
    cs = m_flx.add_ocs_contour_plot(lon,
                                    lat,
                                    ocs_flx[ i, :, : ],
                                    vmin=ocs_flx.min(),
                                    vmax=ocs_flx.max(),
                                    t_str=t_ocs_flx[i],
                                    cmap=cm.get_cmap('Blues_r'),
                                    cbar_t_str='mol m$^{-2}$ s$^{-1}$')
    plt.colorbar(mappable=cs,
                 ax=ax_flx,
                 format='%0.2e')
