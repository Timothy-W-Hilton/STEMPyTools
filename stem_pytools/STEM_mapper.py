import matplotlib.pyplot as plt
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
                 vmax=None):
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
                               fast_or_pretty=fast_or_pretty)
        cm = self.map.add_ocs_contour_plot(self.lon_stem,
                                           self.lat_stem,
                                           self.field_124x124,
                                           n_levs=n_levs,
                                           vmin=vmin,
                                           vmax=vmax)
        plt.colorbar(cm, cax=self.map.ax_cmap)

        return(self)
