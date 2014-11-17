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
    America suitable for plotting STEM 124x124 grid output.
    """

    def __init__(self,
                 t_str="STEM OCS",
                 col_missing=None,
                 use_color=True,
                 map_axis=None,
                 cb_axis=None,
                 missing_axis=None,
                 fig_sz_x=None,
                 fig_sz_y=None):
        """
        class constructor for NAMapFigure.

        INPUT PARAMETERS:
        col_missing: color to use for missing data
        use_color: {True}|False: if True, plot in color.  If false, black and white.
        map_axis: a matplotlib.axes instance.  If provided, the map is
           drawn there.  Default is None, in which case a figure and axes are
           created.
        cb_axis: an axes in which to draw a colorbar (see
           matplotlib.pyplot.colorbar).  May be one of hte following:
              - None (default): no colorbar is drawn
              - a matplotlib.axes instance: the colorbar is drawn here
              - True: an axes is created next to the map axes and the
                colorbar drawn there
        missing_axis: an axes in which the color legend for missing data
           is shown.  Options for specifying are the same as cb_axis.
        fig_sz_x: the horizontal size of the map figure.  Ignored if
           map_axis is specified.
        fig_sz_y: the vertical size of the map figure.  Ignored if
           map_axis is specified.

        Notes:
        If either fig_sz_x or fig_sz_y are unspecified both arguments are
        ignored and the matplotlib default figure size is used.
        """
        if map_axis is None:
            # no axis provided for map, so create a figure and map axis
            if (fig_sz_x is None or fig_sz_y is None):
                fig_sz = None
            else:
                fig_sz = (fig_sz_x, fig_sz_y)
            self.fig = plt.figure(figsize=fig_sz)

            if (cb_axis is None) and (missing_axis is None):
                #no colorbar, so use more of the horizontal extent
                figdim = [0.03, 0.03, 0.94 , 0.94]
            else:
                #save room for colorbar
                figdim = [0.10, 0.05, 0.7 , 0.9]

            self.ax_map = self.fig.add_axes(figdim, frame_on=True)

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
            #self.col_land = '#1C1C1C'
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
        if True:
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
        map.

        INPUT PARMATERS:
        lons: MxN numpy array of longitudes
        lats: MxN numpy array of latitudes
        data: MxN numpy array of data to be contoured
        t_str: title string for the map axes
        vmax: maximum value for color scale.  Default is data.max()
        vmin: minimum value for color scale.  Default is data.min()
        n_levs: number of contour levels.  Default is 20.
        cmap: color map to use for contours.  e.g. cm.get_cmap('Blues').  See
           http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps.
        extend: one of [ 'neither' | 'both' | 'min' | 'max' ], and is
           passed to matplotlib.pyplot.contourf via the 'extend' keyword.
        cbar_t_str: title string for colorbar
        colorbar_args: dict; additional keyword arguments to be passed
            to matplotlib.pypolt.colorbar
        """

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

        #self.draw_mw_box()
        return(cs)

    def add_scatter_plot(self, lons, lats):
        """Plot dots at the locations specified by (lons, lats)."""
        self.map.scatter(lons, lats, latlon=True)

        return(None)

    def draw_mw_box(self, boxcol='black'):
        """
        draw a box around the midwestern region where OCS fluxes were
        enhanced for the pseudo-data
        """
        lat = np.array((40,44))  #deg N
        lon = np.array((-96, -87)) #deg W
        x, y = self.map(
            (lon.min(), lon.min(), lon.max(), lon.max(), lon.min()),
            (lat.min(), lat.max(), lat.max(), lat.min(), lat.min()))
        self.map.plot(x, y, ax=self.ax_map, color=boxcol)
