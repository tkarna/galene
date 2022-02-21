"""
Read Nemo 4.0 station output files and store them as timeprofile files.
"""
import galene as ga

dataset_id = 'run001'
var_list = ['temp', 'psal', 'slev']

search_pattern = '../{:}/run_201*/station_*.nc'.format(dataset_id)

ga.nemo_reader.concatenate_nemo_station_data(search_pattern, dataset_id, var_list)
