"""Provides a number of functions useful for parsing STEM input and
output netcdf files into Python."""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

import na_map

def parse_STEM_AQOUT( aq_fname=None, t0=None, t1=None ):
    """ Parse a OCS from a STEM AQOUT file.  Assumes that the OCS
    variable in the netcdf file is called CO2_TRACER1."""
    aq_out = Dataset( aq_fname, 'r', format='NETCDF4' )
    # read timestamps to datetime.datetime
    t = np.squeeze( aq_out.variables[ 'TFLAG' ] )
    t_dt = [ datetime.strptime( str( this[0] ) + \
                         str(this[1]).zfill(6), '%Y%j%H%M%S' ) \
                         for this in t ]
    # find the requested timestamps
    t_idx = ( t_dt >= t0 ) and ( t_dt <= t1 )
    # retrieve the requested [OCS] data
    ocs = aq_out.variables[ 'CO2_TRACER1' ][ t_idx, 0, :, : ]

def parse_STEM_coordinates( topo_fname ):
    """Parse STEM grid latitude and longitude."""
    topo = Dataset( topo_fname, 'r', format='NETCDF4' )

    lat=np.squeeze( topo.variables['LAT'] )
    lon=np.squeeze( topo.variables['LON'] )

    return( lon, lat )

def initialize_STEM_map2( ax=None ):

    # if no axis provided, create one
    if ax == None:
        fig = plt.figure()
        ax = plt.add_subplot( 1, 1, 1 )

    # define plot corners in projected coordinates.  This pretty much
    # zooms the default view from the Basemap ortho projection so that
    # North America fills the plotting area.  ul is upper left, lr is
    # lower right, etc.  These corners were determined by manually
    # zooming the map produced below.  --TWH
    [ xmin, xmax ] = [ 3.16e6, 9.88e6 ]
    [ ymin, ymax ] = [ 3.07e6, 9.89e6 ]

    # initialize an orthographic projection centered on North America.
    map = Basemap(projection='ortho', lat_0=45, lon_0=-100, resolution='l', \
                  llcrnrx=xmin, llcrnry=ymin, urcrnrx=xmax, urcrnry=ymax )

    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    ## map.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    ## map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.

    # map.drawmeridians(np.arange(0,360,30), labels=[1,0,0,0])
    map.drawparallels(np.arange(-90,90,30), labels=[0,0,0,1 ])
    map.drawmeridians(np.arange(0,360,30) )
    map.drawparallels(np.arange(15,90,20) )

    # set the axes limits to (roughly) North America
    ax = plt.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    return( map, ax )

def initialize_STEM_map():

    return( na_map.NAMapFigure( missing_axis=False ) )

    
def plot_grid_cell_centers( topo_fname, m, ax ):

    lon, lat = parse_STEM_coordinates( topo_fname )

    # plot the grid cell centers
    m.scatter( lon, lat, latlon=True )

    return( ax )
