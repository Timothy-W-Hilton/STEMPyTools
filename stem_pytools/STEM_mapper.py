import matplotlib.pyplot as plt
from mpl_toolkits.basemap import maskoceans

import domain
from stem_pytools.na_map import NAMapFigure


class Mapper124x124(object):
    """Draws a quick map of a 124x124 STEM field superimposed over North
    America.  The 124x124 field is assumed to inhabit the
    stereographic projection grid typically used for N American-wide
    STEM runs in the Campbell lab.

    """

    def __init__(self, field_124x124=None):
        """Class constructor

        INPUTS
        field_124x124: a 124 by 124 numpy array
        """
        d = domain.STEM_Domain()
        self.lat_stem = d.get_lat()
        self.lon_stem = d.get_lon()

        self.field_124x124 = field_124x124

    def draw_map(self,
                 fast_or_pretty='fast',
                 t_str='124x124 quick plot',
                 n_levs=9,
                 vmin=None,
                 vmax=None,
                 cmap=None,
                 norm=None,
                 maskoceans_switch=False,
                 label_latlon=True,
                 cbar_fmt_str=None):
        """draw a quick map of a 124x124 field on the STEM domain

        INPUTS
           fast_or_pretty={'fast'}|'pretty: if 'fast', uses the default
               basemap projection, which renders quickly.  If 'pretty'
               uses the 'satellite' projection, which looks nicer but
               renders more slowly.
           t_str={'124x124 quick plot'}: title for the map
           n_levs: integer {9}; number of contour levels for the map
           vmin: minumum value for color scale; if unspecified uses
              field_124x124.min()
           vmax=: maximum value for color scale; if unspecified uses
              field_124x124.max()
        """

        self.map = NAMapFigure(t_str=t_str,
                               cb_axis=True,
                               fast_or_pretty=fast_or_pretty,
                               label_latlon=label_latlon)

        if maskoceans_switch:
            self.field_124x124 = maskoceans(
                self.lon_stem, self.lat_stem, self.field_124x124,
                resolution='f')

        cm = self.map.map.pcolor(self.lon_stem,
                                 self.lat_stem,
                                 self.field_124x124,
                                 vmin=vmin,
                                 vmax=vmax,
                                 cmap=cmap,
                                 norm=norm,
                                 latlon='True')
        plt.colorbar(cm, cax=self.map.ax_cmap, format=cbar_fmt_str)

        return(self)
