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
    parser.add_argument('-f', '--files', metavar='file', nargs='+',
                        help='a netCDF file to open, multiple values accepted')
    parser.add_argument('--remove-mean', action='store_true',
                        help='remove mean from time series prior to plotting')
    parser.add_argument('-s', '--startdate',
                        help='Start date of time series, in format "YYYY-MM-DDTHH" or "YYYY-MM-DD", e.g. "2006-05-01T00".')
    parser.add_argument('-e', '--enddate',
                        help='Start date of time series, in format "YYYY-MM-DDTHH" or "YYYY-MM-DD", e.g. "2006-05-03T00".')
    parser.add_argument('-v', '--variable', help='Variable to read, e.g. "slev"')
    parser.add_argument('-l', '--location_name', help='Location to read, e.g. "Helsinki"')
    parser.add_argument('-d', '--datasets', help='Comma-separated list of dataset_id to read, e.g. "obs,run01"')
    args = parser.parse_args()
    sd = args.startdate
    if sd is not None:
        sd = dateparse(sd)
    ed = args.enddate
    if ed is not None:
        ed = dateparse(ed)
    cube_list = []
    if args.datasets is not None:
        assert args.location_name is not None, 'location_name must be set'
        assert args.variable is not None, 'variable must be set'
        for dataset_id in args.datasets.split(','):
            dset = so.read_dataset(
                dataset_id, 'timeseries', args.variable,
                location_name=args.location_name,
                start_time=sd, end_time=ed, verbose=False
            )
            c_list = list(dset.values())
            cube_list.extend(c_list)
    if args.files is not None:
        c_list = load_cubes(args.files, start_time=sd, end_time=ed)
        cube_list.extend(c_list)
    if len(cube_list) == 0:
        raise Exception('No timeseries data loaded. Either files or datasets must be defined')
    if args.remove_mean:
        for c in cube_list:
            c.data -= c.data.mean()
    make_interactive_plot(cube_list)


if __name__ == '__main__':
    parse_args()
