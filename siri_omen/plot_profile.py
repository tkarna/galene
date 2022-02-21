"""
Vertical profile plotting routines
"""
import numpy
import matplotlib.pyplot as plt
import os
from . import utility


def plot_profile(ax, cube_list, label_attr='dataset_id', xlim=None,
                 plot_style='line', label_alias=None, title=None,
                 legend_kwargs=None, **kwargs):
    """
    Plots vertical profile objects in the given axes.

    :arg cube_list: list of cube objects to plot. Can be any
        iterable container.
    """

    def get_label(cube, attr_name):
        label = None
        if attr_name in cube.attributes:
            label = cube.attributes.get(attr_name)
        if attr_name == 'date':
            label = str(utility.get_cube_datetime(cube, 0))
        return label

    if isinstance(plot_style, str):
        style_list = [plot_style] * len(cube_list)
    else:
        style_list = plot_style

    for c, style in zip(cube_list, style_list):
        label = get_label(c, label_attr)
        if label_alias is not None:
            label = label_alias.get(label, label)
        if isinstance(c.data, numpy.ma.MaskedArray) \
                and numpy.all(c.data.mask):
            # if all data is bad, skip
            continue
        depth_coord = c.coord('depth')
        if not depth_coord.has_bounds():
            depth_coord.guess_bounds()
        if style == 'cell':
            depth = depth_coord.bounds
            vals = numpy.tile(c.data[:, numpy.newaxis], (1, 2))
            assert vals.shape == depth.shape
            ax.plot(vals.ravel(), depth.ravel(), label=label, **kwargs)
        elif style == 'point':
            kw = dict(kwargs)
            kw['linestyle'] = 'none'
            kw['marker'] = '.'
            ax.plot(c.data, depth_coord.points, label=label, **kw)
        else:
            ax.plot(c.data, depth_coord.points, label=label, **kwargs)
        xlabel = '{:} / {:}'.format(
            c.name().replace('_', ' ').capitalize(), c.units)
        ax.set_xlabel(xlabel)
        ylabel = '{:} / {:}'.format(
            depth_coord.name().capitalize(), depth_coord.units)
        ax.set_ylabel(ylabel)
    ax.grid(True)
    if legend_kwargs is None:
        legend_kwargs = {}
    legend_kwargs.setdefault('loc', 'upper left')
    legend_kwargs.setdefault('bbox_to_anchor', (1.02, 1.0))
    ax.legend(**legend_kwargs)
    if xlim is not None:
        if numpy.array(xlim).size > 2:
            cur_xlim = ax.get_xlim()
            target_xlim = None
            for lim in xlim:
                if cur_xlim[0] >= lim[0] and cur_xlim[1] <= lim[1]:
                    target_xlim = lim
                    break
            if target_xlim is not None:
                ax.set_xlim(target_xlim)
        else:
            ax.set_xlim(xlim)
    # set correct y axis orientation (assuming y is depth, positive downwards)
    ylim = ax.get_ylim()
    ax.set_ylim(sorted(ylim)[::-1])
    if title is None:
        loc_names = utility.unique([c.attributes['location_name']
                                    for c in cube_list])
        date = utility.get_cube_datetime(cube_list[0], 0)
        date_str = date.strftime('%Y-%m-%d')
        titles = utility.unique(loc_names)
        title = ', '.join(titles).strip()
        title += ' ' + date_str
    ax.set_title(title)


def save_profile_figure(cube_list, output_dir=None, plot_root_dir=None,
                        **kwargs):
    """
    Makes a default time series plot and saves it to disk.
    """
    fig = plt.figure(figsize=(5, 12))
    ax = fig.add_subplot(111)

    plot_profile(ax, cube_list, **kwargs)

    imgfile = utility.generate_img_filename(cube_list,
                                            output_dir=output_dir,
                                            root_dir=plot_root_dir,
                                            )
    dir, filename = os.path.split(imgfile)
    utility.create_directory(dir)

    print('Saving image {:}'.format(imgfile))
    fig.savefig(imgfile, dpi=200, bbox_inches='tight')
    plt.close(fig)
