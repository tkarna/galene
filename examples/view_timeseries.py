"""
View netCDF time series in an interactive timeseries plot with bokeh.

Plot opens in a browser window.
"""
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import Category10 as palette
import itertools
import numpy
import argparse
import siri_omen as so
import cftime
from dateutil.parser import parse as dateparse
import datetime

color_cycle = itertools.cycle(palette[10])


def make_interactive_plot(cube_list):

    p = figure(x_axis_type='datetime', plot_width=1400, plot_height=700)

    for cube in cube_list:
        time_coord = cube.coord('time')
        time_units = time_coord.units
        time_array = numpy.array(time_coord.points)
        # get cftime calendar-aware datetime
        datetime_list = time_units.num2date(time_array)
        # convert to python datetime
        datetime_list = [
            datetime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second) for d
            in datetime_list
        ]
        values = cube.data.flatten()

        legend = cube.attributes['dataset_id']
        c = next(color_cycle)
        p.line(datetime_list, values, legend_label=legend, color=c)

    ylabel = '{:} [{:}]'.format(cube.standard_name.replace('_', ' '),
                                cube.units)
    p.yaxis.axis_label = ylabel
    p.title.text = cube.attributes['location_name']

    output_file('timeseries.html')
    show(p)


def load_cubes(file_list, start_time=None, end_time=None):
    cube_list = []
    var_list = ['slev', 'temp', 'psal', 'iceextent', 'icevol']
    for f in file_list:
        c = None
        for v in var_list:
            try:
                c = so.load_cube(f, v,
                                 start_time=start_time, end_time=end_time)
                break
            except Exception:
                pass
        assert c is not None, 'Could not open file {:}'.format(f)
        cube_list.append(c)
    return cube_list


def parse_args():

    parser = argparse.ArgumentParser(
        description='Open netCDF time series files in an interactive plot.',
        # includes default values in help entries
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('files', metavar='file', nargs='+',
                        help='a netCDF file to open, multiple values accepted')
    parser.add_argument('--remove-mean', action='store_true',
                        help='remove mean from time series prior to plotting')
    parser.add_argument('-s', '--startdate',
                        help='Start date of time series, in format "YYYY-MM-DDTHH" or "YYYY-MM-DD", e.g. "2006-05-01T00".')
    parser.add_argument('-e', '--enddate',
                        help='Start date of time series, in format "YYYY-MM-DDTHH" or "YYYY-MM-DD", e.g. "2006-05-03T00".')
    args = parser.parse_args()
    sd = dateparse(args.startdate)
    ed = dateparse(args.enddate)
    cube_list = load_cubes(args.files, start_time=sd, end_time=ed)
    if args.remove_mean:
        for c in cube_list:
            c.data -= c.data.mean()
    make_interactive_plot(cube_list)


if __name__ == '__main__':
    parse_args()
