"""
Make 2D map plots.
"""
import numpy
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
from .plot_taylor_diag import get_point_style_cycler
from matplotlib.tri import Triangulation


class GeographicPlot(object):
    """
    A map plot object.
    """
    def __init__(self, ax=None, fig=None, rect=None, projection=None,
                 extent=[6, 32, 53, 66], symmetric_colorbar=False,
                 draw_grid=True, grid_kwargs={}, grid_zorder=2):
        self.symmetric_colorbar = symmetric_colorbar
        self.val_max_magnitude = 0.0

        if projection is None:
            projection = ccrs.Mercator(
                central_longitude=20., min_latitude=0.,
                max_latitude=80., latitude_true_scale=60.0
            )

        if ax is None:
            # create new axes with projection
            self.fig = plt.figure(figsize=(12, 12)) if fig is None else fig
            if rect is None:
                rect = 111

            ax = self.fig.add_subplot(rect, projection=projection)
            self.ax = ax
        else:
            # use user-provided axis, must have a geographic projection
            self.ax = ax
            self.fig = self.ax.figure

        if extent is not None:
            ax.set_extent(extent)

        if draw_grid:
            grid_kwargs.setdefault('linewidth', 0.4)
            gl = ax.gridlines(draw_labels=True, zorder=grid_zorder,
                              **grid_kwargs)
            gl.top_labels = False
            gl.right_labels = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            self.gl = gl

        self.point_style = get_point_style_cycler()
        self.point_style_iter = iter(self.point_style)

        self.sample_points = []

    def add_station(self, x, y, label, label_to_legend=False, textargs=None,
                    transform=None, hide_label=False, **kwargs):
        if transform is None:
            transform = ccrs.PlateCarree()
        extent = self.ax.get_extent(transform)
        if (x < extent[0] or x > extent[1] or y < extent[2] or y > extent[3]):
            print(f'Station {label} {x},{y} outside the plot, skipping.')
            return None, None
        if label_to_legend and not hide_label:
            kwargs.setdefault('label', label)
        else:
            kwargs.pop('label', None)
        kwargs.setdefault('zorder', 3)
        sty = next(self.point_style_iter)
        sty.update(kwargs)
        p = self.ax.plot(x, y, transform=transform, **sty)
        if label_to_legend and not hide_label:
            self.sample_points.append(p)
            txt = None
        elif hide_label:
            txt = None
        else:
            if textargs is None:
                textargs = {}
            textargs.setdefault('horizontalalignment', 'left')
            textargs.setdefault('verticalalignment', 'top')
            textargs.setdefault('fontsize', 8)
            offset = 0.1
            ha = textargs['horizontalalignment']
            if 'ha' in textargs:
                ha = textargs['ha']
            if ha == 'right':
                offset *= -1
            txt = self.add_text(x + offset, y, label,
                                transform=transform, **textargs)
        return p, txt

    def add_text(self, x, y, text, transform=None, **kwargs):
        if transform is None:
            transform = ccrs.PlateCarree()
        txt = self.ax.text(x, y, text, transform=transform, **kwargs)
        return txt

    def add_unstructured_mesh(self, x=None, y=None, connectivity=None,
                              tri=None, transform=None, **kwargs):
        if transform is None:
            transform = ccrs.PlateCarree()
        has_xyc = x is not None and y is not None and connectivity is not None
        has_tri = tri is not None
        msg = 'either x, y, connectivity or tri arguments must be provided'
        assert (has_tri and not has_xyc) or (has_xyc and not has_tri), msg

        kwargs.setdefault('linewidth', 0.3)
        kwargs.setdefault('zorder', 2)

        if not has_tri:
            tri = Triangulation(x, y, connectivity)

        p = self.ax.triplot(tri, transform=transform, **kwargs)
        return p

    def add_unstructured_field(self, x=None, y=None, connectivity=None,
                               tri=None, values=None, transform=None,
                               levels=31,
                               **kwargs):
        if transform is None:
            transform = ccrs.PlateCarree()
        has_xyc = x is not None and y is not None and connectivity is not None
        has_tri = tri is not None
        msg = 'either x, y, connectivity or tri arguments must be provided'
        assert (has_tri and not has_xyc) or (has_xyc and not has_tri), msg
        assert values is not None, 'values must be provided'

        kwargs.setdefault('zorder', 2)

        if not has_tri:
            tri = Triangulation(x, y, connectivity)

        p = self.ax.tricontourf(tri, values, levels, transform=transform,
                                **kwargs)
        return p

    def add_cube(self, cube, kind='pcolormesh', **kwargs):
        if self.symmetric_colorbar:
            max_mag = max(abs(cube.data.min()), abs(cube.data.max()))
            self.val_max_magnitude = max(max_mag, self.val_max_magnitude)
            vmin = -self.val_max_magnitude
            vmax = self.val_max_magnitude
            kwargs.setdefault('vmin', vmin)
            kwargs.setdefault('vmax', vmax)

        kwargs.setdefault('zorder', 2)
        if kind == 'pcolormesh':
            p = iplt.pcolormesh(cube, axes=self.ax, **kwargs)
        elif kind == 'pcolor':
            p = iplt.pcolor(cube, axes=self.ax, **kwargs)
        elif kind == 'contourf':
            p = iplt.contourf(cube, axes=self.ax, **kwargs)
        elif kind == 'contour':
            p = iplt.contour(cube, axes=self.ax, **kwargs)
        else:
            raise RuntimeError('Unknown plot kind: {:}'.format(kind))
        return p

    def add_quiver(self, cube_x, cube_y, c=None, transform=None, **kwargs):
        x = cube_x.coord('longitude').points
        y = cube_x.coord('latitude').points
        u = cube_x.data
        v = cube_y.data
        if transform is None:
            transform = ccrs.PlateCarree()
        assert u.shape == v.shape
        kwargs.setdefault('zorder', 3)
        kwargs.setdefault('pivot', 'tail')
        kwargs.setdefault('headwidth', 4.5)
        kwargs.setdefault('headlength', 6)
        if c is not None:
            p = self.ax.quiver(x, y, u, v, c, transform=transform, **kwargs)
        else:
            p = self.ax.quiver(x, y, u, v, transform=transform, **kwargs)
        return p

    def add_feature(self, feature_name, scale='50m', **kwargs):
        kwargs.setdefault('zorder', 1)
        if feature_name == 'land':
            facecolor = kwargs.pop('facecolor', '0.8')
            edgecolor = kwargs.pop('edgecolor', 'none')
            land_feature = cfeature.NaturalEarthFeature(
                'physical', 'land', scale,
                edgecolor=edgecolor, facecolor=facecolor)
            self.ax.add_feature(land_feature, **kwargs)
        elif feature_name == 'coastlines':
            self.ax.coastlines(resolution='50m', **kwargs)

    def add_legend(self, **kwargs):
        """
        Add a legend to the plot
        """
        kwargs.setdefault('prop', dict(size='small'))
        kwargs.setdefault('loc', 'upper left')
        kwargs.setdefault('bbox_to_anchor', (1.01, 1.0))
        nsamples = len(self.sample_points)
        ncolumns = int(numpy.ceil(float(nsamples) / 40))
        kwargs.setdefault('ncol', ncolumns)
        self.ax.legend(numpoints=1, **kwargs)

    def add_colorbar(self, p, label=None, location='left', width=0.02, pad=0.03,
                     cax=None, **kwargs):
        if cax is None:
            # create new axes for colorbar
            orientation = 'vertical' if location in [
                'right', 'left'] else 'horizontal'
            kwargs.setdefault('orientation', orientation)
            pos = self.ax.get_position().bounds
            if location == 'left':
                x = pos[0] + pos[2] + pad * pos[2]
                cax = self.fig.add_axes([x, pos[1], width, pos[3]])
            elif location == 'bottom':
                y = pos[1] - pad * pos[3] - width
                cax = self.fig.add_axes([pos[0], y, pos[2], width])
        cb = plt.colorbar(p, cax=cax, **kwargs)
        if label is not None:
            cb.set_label(label)
        return cb

    def add_title(self, title):
        """
        Add a title to the plot
        """
        self.ax.set_title(title)
