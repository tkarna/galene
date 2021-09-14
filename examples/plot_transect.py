import numpy
import matplotlib.pyplot as plt
import siri_omen as so
from siri_omen.numpy_interface import create_transect_cube
import iris

npoints = 40
ndepth = 12

# generate coordinates
lat = numpy.linspace(61.0, 61.5, npoints)
lon = numpy.linspace(20.5, 19.5, npoints)
depth = numpy.linspace(0, 40, ndepth)
time = numpy.linspace(0, 5*24*3600., npoints)

# generate artificial data
distance = numpy.linspace(0, 1, npoints)
a, b = numpy.meshgrid(depth, distance, indexing='ij')
values = numpy.sin(10*b) * a/20

# add required metadata
attributes = {}
attributes['dataset_id'] = 'testdata'  # e.g. mission name, uivelo-20211
attributes['location_name'] = 'transectA'  # e.g. transect name, segment4-5

cube = create_transect_cube(lat, lon, depth, values, time,
    'sea_water_practical_salinity', '1e-3', attributes=attributes)
print(cube)

# save in standard format
so.save_cube(cube)

# plot
so.save_timetransect_figure([cube], None)
plt.show()
