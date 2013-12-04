import os.path
import numpy as np
import numpy.ma as ma
from datetime import datetime
from mpl_toolkits.basemap import maskoceans

import STEM_vis

topo_fname = '/Users/tim/work/Data/STEM/input/TOPO-124x124.nc'
aq_fname = os.path.join('/Users/tim/work/Data/STEM/output/',
                        'AQOUT-124x124-22levs-casa-cos_2008_2009.nc')

[ ocs, t ] = STEM_vis.parse_STEM_AQOUT(aq_fname,
                                       datetime(2008, 6, 1),
                                       datetime(2008, 6, 8))
lon, lat = STEM_vis.parse_STEM_coordinates( topo_fname )

mask = maskoceans( lon, lat, np.zeros( lon.shape ), inlands = False )
mask = mask.mask  # get rid of the dummy zeros
mask = np.tile( mask, [ ocs.shape[0], 1, 1 ] )
ocs = ma.masked_array( ocs, mask = mask )

t_lab = None #initialize time label
#for i in range( ocs.shape[0] ):
for i in range( 13,14 ):
    t0 = datetime.now()
    m = STEM_vis.initialize_STEM_map()
    cs = m.add_ocs_contour_plot(lon,
                                lat,
                                ocs[ i, :, : ],
                                t_str=t[i])

    fname = datetime.strftime( t[i], './Maps3/%Y-%m-%dT%H%M_STEM_OCS.pdf' )
    # print 'saving ' + fname + '(' + str(datetime.now() - t0) + ')'
    # plt.savefig( fname, bbox_inches = 0 )
    # plt.close( m.fig )

## TODO
##    - documentation
