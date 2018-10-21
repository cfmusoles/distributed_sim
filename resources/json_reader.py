import json
from pprint import pprint

json_file = 'multiareaModel.json'#'viscortex_raw_data.json' 'viscortex_processed_data.json' 
key_attribute = 'synapses'#'surface_data'  'SLN_completed'

# Model.py has the algorithm to build the model ('default_Data_Model_.json')
# Copy viscortex_raw_data.json. Model.py and VisualCortex_Data.py to the parent directory
# Remove __init__.py file from multiarea_model folder
# run python Model.py
# May have to remove from config import base_path and create a dummy base_path=''
#population_list --> list of pop references (23E, 23I...)
#area_list --> list of area references
#synapse_weights_mean: E 87.8085, I -1404.94 (E * -16 = I)
#synapse_weights_sd: 0.1 of the mean
#neuron_numbers --> per population, number of neurons per layer
#synapses == realistic_synapses --> number of synapses from to any layers (populations)
    # [target_area][target_pop][source_area][source_pop]


with open(json_file) as f:
    data = json.load(f)

for key in data.keys():
    pprint(key)

pprint(data[key_attribute])