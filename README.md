# galene

Post-processing tools for ocean model outputs.


## Installation

### Install Iris

`galene` uses [Iris](https://scitools-iris.readthedocs.io/en/stable/) as it's
data model. Iris is easiest to install with Anaconda Python. These instructions have been tested for Anaconda version `Anaconda3-2019.03`.

It recommended to create a specific conda environment for Iris, for example:
```bash
conda create --name iris python=3
```

The environment must be activated prior usage:
```bash
conda activate iris
```

Then we can install Iris. See the Iris website for up-to-date [installation instructions](https://scitools-iris.readthedocs.io/en/stable/installing.html).

```bash
conda install -c conda-forge iris
```

### Install galene with pip

Once Iris is installed, and the `iris` environment is active, install
`galene` with

```bash
pip install -e /path/to/galene/
```

## Features

- time series plots
- vertical profile plots
- comparison of datasets, interpolation on common grid
- computation of statistics
- Taylor and target diagrams
- Geographical plots (maps)

## Data model

`galene` data model uses Iris Cube objects to represent data
(see [Iris documentation](https://scitools-iris.readthedocs.io/en/stable/userguide/iris_cubes.html)).
In addition, two metadata entries are required:

1. `cube.attributes['dataset_id']`: A string that identifies the data set,
    e.g. `mynemorun1` or `observations`
2. `cube.attributes['location_name']`: A string that identifies the spatial
    location of the data, e.g. station or transect identifier.

Currently these kinds of geospatial data are supported (dimension coordinates in parentheses):

- timeseries: (time) + auxiliary scalar depth
- surfacetrack: (time) + auxiliary latitude, longitude, scalar depth
- profile: (depth) + auxiliary scalar latitude, longitude
- timeprofile: (time, depth)
- transect: (depth, index) + auxiliary latitude, longitude, depth
- timetransect: (time, depth, index)

All data types are stored as a `Cube` objects. The function `get_cube_datatype(cube)` returns the data type as a string.

## Reading netCDF data

`galene` can read netCDF files that contain sufficient metadata.
It is recommended to first generate Iris Cube objects and then store them to
disk.

To support reading in various model output files, some metadata editing may be
necessary. See `galene/nemo_reader.py`.
