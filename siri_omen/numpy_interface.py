import iris
import cf_units
from collections import Iterable


def create_timeseries_cube(time, values, standard_name, units,
                           depth=0.0,
                           time_units='seconds since 1970-01-01 00:00:00-00',
                           calendar='gregorian',
                           attributes=None):
    """
    Create a timeseries Cube object.

    Parameters
    ----------

    time: array_like (ntime, )
        array of time stamps
    values: array_like (ntime, )
        array of values
    standard_name: str
        CF standard name of the variable.
        E.g. 'sea_water_temperature'
    units: srt
        CF standard units of the variable.
        E.g. 'degree_C'
    depth: float
        Depth in meters. Default: 0.0
    time_units: str
        CF standard units of the time variable.
        E.g. 'seconds since 1970-01-01 00:00:00-00'
    calendar: str
        CF standard calendar of the time variable.
        Default: 'gregorian'
    attributes: dict
        Dictionary of metadata attributes.

    Returns
    -------

    Timeseries as an Iris Cube object.

    """
    msg = 'length of time and values arrays does not match.'
    assert len(time) == len(values), msg

    # define time coordinate
    time_units_obj = cf_units.Unit(time_units, calendar=calendar)

    # make time DimCoord
    time_coord = iris.coords.DimCoord(time, standard_name='time',
                                      units=time_units_obj)

    # make a cube
    cube = iris.cube.Cube(values, standard_name=standard_name, units=units)

    # add time coordinate to dimension 0
    cube.add_dim_coord(time_coord, 0)

    # add aux depth coordinate (required for time series)
    dep_coord = iris.coords.DimCoord([depth], standard_name='depth',
                                     units='m')
    cube.add_aux_coord(dep_coord)

    if attributes is not None:
        # add attributes
        cube.attributes.update(attributes)

    return cube


def create_transect_cube(latitude, longitude, depth, values, time,
                         standard_name, units,
                         time_units='seconds since 1970-01-01 00:00:00-00',
                         calendar='gregorian',
                         attributes=None):
    """
    Create a transect Cube object. Transect is a 2D data set with dimensions
    (depth, points) where points stand for points in the horizontal plane.

    Parameters
    ----------

    latitude, longitude: array_like (npoints, )
        transect latitude, longitude coordinates
    depth: array_like (ndepth, )
        transect depth values
    values: array_like (ndepth, npoints)
        array of values
    time: scalar or array_like (npoints, )
        transect time stamp(s)
    standard_name: str
        CF standard name of the variable.
        E.g. 'sea_water_temperature'
    units: srt
        CF standard units of the variable.
        E.g. 'degree_C'
    depth: float
        Depth in meters. Default: 0.0
    time_units: str
        CF standard units of the time variable.
        E.g. 'seconds since 1970-01-01 00:00:00-00'
    calendar: str
        CF standard calendar of the time variable.
        Default: 'gregorian'
    attributes: dict
        Dictionary of metadata attributes.

    Returns
    -------

    Transect as an Iris Cube object.

    """
    assert len(values.shape) == 2, 'values must be a 2D array'
    ndepth, npoints = values.shape
    assert len(latitude) == npoints, \
        'latitude length must match value array shape[1]'
    assert len(longitude) == npoints, \
        'longitude length must match value array shape[1]'
    assert len(depth) == ndepth, 'depth length must match value array shape[0]'

    # make a cube
    cube = iris.cube.Cube(values, standard_name=standard_name, units=units)

    # define time coordinate
    time_units_obj = cf_units.Unit(time_units, calendar=calendar)
    time_coord = iris.coords.AuxCoord(time, standard_name='time',
                                      units=time_units_obj)

    # add time coordinate
    if isinstance(time, Iterable):
        cube.add_aux_coord(time_coord, 1)
    else:
        cube.add_aux_coord(time_coord)

    # add aux coords
    dep_coord = iris.coords.AuxCoord(depth, standard_name='depth', units='m')
    cube.add_aux_coord(dep_coord, 0)

    lat_coord = iris.coords.AuxCoord(latitude, standard_name='latitude',
                                     units='degree_north')
    cube.add_aux_coord(lat_coord, 1)
    lon_coord = iris.coords.AuxCoord(longitude, standard_name='longitude',
                                     units='degree_east')
    cube.add_aux_coord(lon_coord, 1)

    if attributes is not None:
        # add attributes
        cube.attributes.update(attributes)

    return cube
