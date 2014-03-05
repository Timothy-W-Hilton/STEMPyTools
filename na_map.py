""" Module to provide a Matplotlib figure containing a map of North
America suitable for plotting STEM 124x124 grid output"""

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import pdb

R_EARTH = 6371007.181000  #Earth radius in meters

class NAMapFigure(object):
    """ Class to provide a Matplotlib figure containing a map of North
    America suitable for plotting STEM 124x124 grid output"""

    def __init__(self,
                 t_str="STEM OCS",
                 col_missing=None,
                 use_color=True,
                 map_axis=None,
                 cb_axis=None,
                 missing_axis=None,
                 fig_sz_x=None,
                 fig_sz_y=None):

        if map_axis is None:
            # no axis provided for map, so create a figure and map axis
            if (fig_sz_x is None or fig_sz_y is None):
                fig_sz = None
            else:
                fig_sz = (fig_sz_x, fig_sz_y)
            self.fig = plt.figure(figsize=fig_sz)

            self.ax_map = self.fig.add_axes([0.10, 0.05, 0.7 , 0.9],
                                            frame_on=True)

            self.ax_cmap = cb_axis
            if ((cb_axis is not None) and
                (type(self.ax_cmap) is not mpl.axes.Axes)):
                self.ax_cmap = self.fig.add_axes([0.85, 0.15, 0.05, 0.70],
                                                 frame_on=True )
            self.ax_miss = missing_axis
            if ((missing_axis is not None) and
                (type(self.ax_miss) is not mpl.axes.Axes)):
                self.ax_miss = self.fig.add_axes([0.81, 0.1, 0.05, 0.05],
                                                 frame_on=False)
        else:
            # a map axis was provided -- draw the map there
            self.ax_map = map_axis
            self.ax_cmap = cb_axis
            self.ax_miss = missing_axis
            self.fig = self.ax_map.figure

        if t_str is not None:
            self.ax_map.set_title(t_str)


        mapwidth = 8.0e6  # not sure of units for width/height
        mapheight = 6.5e6
        self.map = Basemap(width=mapwidth,
                           height=mapheight,
                           projection='aeqd',
                           lat_0=54,
                           lon_0=-105,
                           resolution='l',
                           area_thresh=1000, #show features larger than 1000 km
                           rsphere=R_EARTH,
                           ax=self.ax_map,
                           fix_aspect=True)

        #define some colors
        if use_color:
            self.col_water = '#B9D3EE'  #SlateGray2
            self.col_land = '#FFF8DC'  #cornsilk
            self.col_land = '#1C1C1C'
            self.col_states = "#0A0A0A"
            self.map_grid_col = "#000000"  #color for map grid lines
        else:
            self.col_water = "#FFFFFF"
            self.col_land = '#EEEEEE'
            self.col_states = "#CCCCCC"
            self.map_grid_col = "#000000"
        if col_missing is None:
            # from ColorBrewer
            self.col_missing = np.array([217.0, 95.0, 2.0, 0.0]) / 255.0
        else:
            self.col_missing = col_missing

        self.map.drawmapboundary(fill_color = self.col_water)
        self.map.drawcoastlines(linewidth=0.5)
        if False:
            self.map.drawstates(color=self.col_states)
            self.map.drawcountries(color=self.col_states)
        # zorder very important - continent fill needs to be at the
        # back!  Otherwise it is drawn on top of everything else
        self.map.fillcontinents(color=self.col_land,
                                lake_color=self.col_water,
                                zorder=0)
        self.map.drawmeridians( meridians=range(0, -180, -20),
                                labels=(0, 0, 0, 1),  # labels on bottom
                                color=self.map_grid_col)
        self.map.drawparallels(circles=range(0, 90, 20),
                               labels=(1, 0, 0, 0),  # labels on left
                               color=self.map_grid_col)

    def add_ocs_contour_plot(self,
                             lons,
                             lats,
                             data,
                             t_str=None,
                             vmax=None,
                             vmin=None,
                             n_levs=20,
                             cmap=cm.get_cmap('Blues'),
                             extend='neither',
                             cbar_t_str=None,
                             colorbar_args={}):
        """Draw filled contours of the specified OCS data over the
        map.  Extend must be one of [ 'neither' | 'both' | 'min' |
        'max' ], and is passed to matplotlib.pyplot.contourf via the
        'extend' keyword."""

        if vmin is None:
            vmin = data.min()
        if vmax is None:
            vmax = data.max()

        contour_levs = np.linspace(vmin,
                                   vmax,
                                   n_levs,
                                   endpoint=True)

        cs = self.map.contourf(lons,
                               lats,
                               data,
                               ax=self.ax_map,
                               levels=contour_levs,
                               extend=extend,
                               latlon=True,
                               cmap=cmap)
        if self.ax_cmap is not None:
            # plot a color legend
            plt.colorbar(mappable=cs,
                         cax=self.ax_cmap,
                         **colorbar_args)
            self.ax_cmap.set_title( cbar_t_str )

        if t_str is not None:
            # place a time label at (140, 20N) (out in the Pacific
            # Ocean west of Mexico)
            t_lab_lon = -140
            t_lab_lat = 20
            #t_str = datetime.strftime( t_str, '%d %B %Y %H:%M' )
            self.ax_map.text(*self.map( t_lab_lon, t_lab_lat ),
                             s=t_str,
                             bbox=dict(facecolor='white', alpha=0.5))
        return(cs)

    def add_scatter_plot(self, lons, lats):
        """Plot dots at the locations specified by (lons, lats)."""
        self.map.scatter(lons, lats, latlon=True)

        return(None)
