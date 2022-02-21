import iris
import galene as ga
import datetime
from scipy.spatial import cKDTree
import numpy
import more_itertools
import glob

var_to_nemo_var = {
    'temp': 'sea_surface_temperature',
    'psal': 'sea_surface_salinity',
}


class NearestNeighborLatLon:
    """
    Find nearest neighbors in (lat,lon) grid taking into account the land mask.
    """
    def __init__(self, cube):
        """
        Generate a nearest neighbor finder from a 2D NEMO field.
        """
        self.lat = cube.coord('latitude').points
        self.lon = cube.coord('longitude').points
        self.landmask = cube.data.mask
        assert len(self.landmask.shape) == 2
        self.lon_2d, self.lat_2d = numpy.meshgrid(self.lon, self.lat,
                                                  copy=True)
        self.wetmask = numpy.nonzero(~self.landmask.ravel())[0]
        # 1D arrays of all wet points
        self.valid_lon = self.lon_2d.ravel()[self.wetmask]
        self.valid_lat = self.lat_2d.ravel()[self.wetmask]
        # make nn tree
        coords = numpy.vstack((self.valid_lon, self.valid_lat)).T
        self.tree = cKDTree(coords)

    def search(self, lon, lat):
        """
        Find the nearest (i,j) index in the 2D grid.

        (i,j) correspond to the first 2 spatial indices of the NEMO grid, i.e.
        latitude and longitude index, respectively.
        """
        xy = numpy.vstack([lon, lat]).T
        dist, index = self.tree.query(xy, k=1)
        index = self.wetmask[index]
        i, j = numpy.unravel_index(index, self.landmask.shape)
        return i, j


def extract_track_from_cube(nemo_cube, track_cube, time_pad, dataset_id,
                            nn_finder=None):
    """
    Extract surface track from NEMO 2d cube
    """
    # crop track time
    st = ga.get_cube_datetime(nemo_cube, 0)
    et = ga.get_cube_datetime(nemo_cube, -1)
    # NOTE do not include start instant to have non-overlapping windows
    target = ga.constrain_cube_time(
        track_cube, st - time_pad, et + time_pad, include_start=False
    )

    def find_nearest_index(src, dst, coord_name):
        src_arr = src.coord(coord_name).points
        dst_arr = dst.coord(coord_name).points
        time_tree = cKDTree(src_arr[:, numpy.newaxis])
        d, index = time_tree.query(dst_arr[:, numpy.newaxis], k=1)
        return index

    if nn_finder is None:
        nn_finder = NearestNeighborLatLon(nemo_cube[0, ...])

    target_lat = target.coord('latitude').points
    target_lon = target.coord('longitude').points

    i_lat, i_lon = nn_finder.search(target_lon, target_lat)
    ntime = len(nemo_cube.coord('time').points)
    if ntime == 1:
        i_time = numpy.zeros_like(i_lat)
    else:
        i_time = find_nearest_index(nemo_cube, target, 'time')

    values = nemo_cube.data[i_time, i_lat, i_lon]

    sname = ga.nemo_reader.map_nemo_sname_to_standard[nemo_cube.standard_name]
    cube = iris.cube.Cube(values, standard_name=sname, units=nemo_cube.units)

    # copy coordinates
    cube.add_dim_coord(target.coord('time'), 0)
    cube.add_aux_coord(target.coord('latitude'), 0)
    cube.add_aux_coord(target.coord('longitude'), 0)
    cube.add_aux_coord(target.coord('depth'))
    for coord_name in ['time', 'latitude', 'longitude', 'depth']:
        cube.coord(coord_name).attributes = {}  # discard coord attributes
    # add attributes
    cube.attributes['location_name'] = target.attributes['location_name']
    cube.attributes['dataset_id'] = dataset_id

    return cube


def extract_from_single_nemo_file(src_file, track_cube, var, dataset_id,
                                  chunksize=50):
    """
    Extract tract from a NEMO 2D output file.
    """
    nemo_var = var_to_nemo_var[var]
    print(f'Reading {nemo_var} from {src_file}')
    nemo_cube = ga.nemo_reader.load_nemo_output(src_file, nemo_var)
    # compute time window for the model data
    # the time range is
    # [first_time_stamp - time_step/2, last_time_stamp + time_step/2]
    # where the latter bound is inclusive
    time_units = nemo_cube.coord('time').units
    time_array = nemo_cube.coord('time').points
    time_step = time_array[1] - time_array[0]
    assert 'seconds' in str(time_units), 'model time should be in seconds'
    time_pad = datetime.timedelta(seconds=numpy.round(time_step/2))

    track_cube.coord('time').convert_units(time_units)

    nn_finder = NearestNeighborLatLon(nemo_cube[0, ...])

    cube_list = []
    for chunk in more_itertools.chunked(range(len(time_array)), chunksize):
        try:
            src = nemo_cube[chunk, ...]
            c = extract_track_from_cube(
                src, track_cube, time_pad, dataset_id, nn_finder=nn_finder
            )
            cube_list.append(c)
        except AssertionError as e:
            pass
    if len(cube_list) == 0:
        return None
    cube = ga.concatenate_cubes(cube_list)
    # sanity checks, all values should be valid
    assert not cube.data.mask.any(), 'Some extracted values are masked'
    assert numpy.isfinite(cube.data).all(), 'Some extracted values are invalid'
    assert numpy.abs(cube.data).max() < 1e10, \
        'Some extracted values are too large'
    return cube


def extract_track(nemo_file_pattern, track_file, var_list, dataset_id, chunksize=50):
    """
    Extract surface track from a series of NEMO output files.
    """
    track_cube = ga.load_cube(dst_file)

    for var in var_list:
        file_list = glob.glob(nemo_file_pattern)
        assert len(file_list) > 0, f'No files found: {nemo_file_pattern}'
        cube_list = []
        for f in sorted(file_list):
            c = extract_from_single_nemo_file(
                f, track_cube, var, dataset_id, chunksize=chunksize
            )
            if c is not None:
                cube_list.append(c)
        cube = ga.concatenate_cubes(cube_list)
        ga.save_cube(cube)


# variables to extract
var_list = ['temp', 'psal']
# name of the run
dataset_id = 'run31b'
# pattern where to search for 2D NEMO output files
src_file_pattern = '../run31b/run_*/output/NORDIC_1h_SURF_grid_T_*.nc'
# Observation track which defines the locations and times (iris cube)
dst_file = 'cmems-ferrybox/surfacetrack/TransPaper/psal/strack_TransPaper_d0.00m_cmems-ferrybox_psal_2014-01-16_2016-11-27.nc'
extract_track(src_file_pattern, dst_file, ['psal'], dataset_id, chunksize=50)
