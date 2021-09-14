import iris
import cf_units


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

    Timeseries as a Iris Cube object.

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
