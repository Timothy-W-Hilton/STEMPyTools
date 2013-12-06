import matplotlib as mpl
mpl.use('Agg')

import os.path
import numpy as np
import numpy.ma as ma
from datetime import datetime
from mpl_toolkits.basemap import maskoceans
import matplotlib.pyplot as plt

import STEM_vis

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
for i in range( 39,40 ):
    t0 = datetime.now()
    m = STEM_vis.initialize_STEM_map()
    cs = m.add_ocs_contour_plot(lon,
                                lat,
                                ocs[ i, :, : ],
                                t_str=t[i])

    # fname = datetime.strftime( t[i], '/home/thilton/Plots/STEM_Maps/%Y-%m-%dT%H%M_STEM_OCS.png' )
    # print 'saving ' + fname + '(' + str(datetime.now() - t0) + ')'
    # plt.savefig( fname )
    # plt.close( m.fig )

## TODO
##    - documentation

    m_flx = STEM_vis.initialize_STEM_map()
    cs = m.add_ocs_contour_plot(lon,
                                lat,
                                ocs_flx[ i, :, : ],
                                t_str=t_ocs_flx[i])
