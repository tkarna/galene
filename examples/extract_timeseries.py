"""
Read NEMO timeprofiles and extract time series at observation depths.
"""
import galene as ga
import iris

var_list = ['temp', 'psal']
obs_id = 'cmems-nrt'
dataset_id = 'run001'

for var in var_list:
    dataset = ga.read_dataset(dataset_id, 'timeprofile', var)
    obs_dataset = ga.read_dataset(obs_id, 'timeseries', var)

    pairs = ga.find_station_pairs(obs_dataset, dataset)

    for key in pairs:
        o = pairs[key][obs_id]
        m = pairs[key][dataset_id]

        def contains_coord(cube, standard_name):
            for c in cube.coords():
                if c.standard_name == standard_name:
                    return True
            return False

        depth = o.coord('depth').points[0]
        if contains_coord(m, 'depth'):
            # interpolate depth
            target = [('depth', depth)]
            cube = m.interpolate(target, iris.analysis.Linear())
        else:
            # assume already at correct depth
            cube = m
            dep_dim = iris.coords.DimCoord(depth,
                                           standard_name='depth',
                                           units='m')
            cube.add_aux_coord(dep_dim, None)
        # drop lon,lat dims
        cube = ga.drop_singleton_dims(cube)

        ga.save_cube(cube)
