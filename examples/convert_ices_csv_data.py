import galene as ga

input_files = ['1231022b.csv', ]
dataset_id = 'ices-ctd'
station_file = 'helcom_stations.csv'

for f in input_files:
    ga.ices_reader.process_cvs_file(f, station_file, dataset_id)
