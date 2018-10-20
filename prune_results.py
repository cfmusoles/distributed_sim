### Python helper file to select experiment dataset (reduction) from files
### The reduction can be done based on any column value from test results file
### based on a comparator (max or min, edit line 55)
### Only experiments with the same seed are compared and reduced
### if there are two groups of results with a seed each, the resulting pruned files will contain 2 results only

import numpy as np
import sys
import os

if len(sys.argv) < 6:
	print("Input error: usage must be python prune_results.py filename base_column min_process geometric_step num_experiments")
	exit()

seed_column = 10 # which column includes seed info on test results
test_name = sys.argv[1]
base_column = int(sys.argv[2])
min_processes = int(sys.argv[3])
geometric_step = int(sys.argv[4])
num_experiments = int(sys.argv[5])

def get_data_from_csv_irregular_columns(filename,delimiter=",",skip_header=True):
    import csv
    datafile = open(filename, 'r')
    datareader = csv.reader(datafile)
    if(skip_header):
        next(datareader,None)
    data = []
    for row in datareader:
        data.append( [ float(elem) for elem in delimiter.join(row).split(delimiter) ] )
    return data

def select_rows_from_csv_irregular_columns(selected_rows,input_filename,output_filename,with_header=True):
    import csv
    datafile = open(input_filename, 'r')
    writebuffer = open(output_filename,'w')
    if with_header:
        writebuffer.write(datafile.readline())
    for counter, row in enumerate(datafile):
        if counter in selected_rows:
            writebuffer.write(row)

# Open test results file
# initialise selected rows as []
# For each unique seed number (numpy.unique)
#   get all rows with that seed number (numpy.in1d)
#   compare the values of base column to select the candidate row (argmin())
#   add row to selected rows
# call select_rows_from_csv on all the file types, passing the selected rows as argument

process_counts = [min_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]

for num_processes in process_counts:
    test_results_filename = test_name + "__" + str(num_processes)
    test_results_filename_output = test_name + "_pruned__" + str(num_processes)
    process_communication_filename = test_name + "_" + str(num_processes)
    trace_filename = test_name + "_" + str(num_processes)

    selected_rows = []
    test_results_data = np.array(get_data_from_csv_irregular_columns(test_results_filename))
    all_seeds = np.unique(test_results_data[:,seed_column])
    for seed in all_seeds:
        selected_data = [[counter,candidate] for counter,candidate in enumerate(test_results_data) if candidate[seed_column] == seed]
        indices = np.array(list(zip(*selected_data))[0])
        candidates = np.array(list(zip(*selected_data))[1])[:,base_column]
        sel_row = candidates.argmin()
        selected_rows.append(indices[sel_row])
        

    # test results file selection
    if os.path.isfile(test_results_filename):
        select_rows_from_csv_irregular_columns(selected_rows,test_results_filename,test_results_filename_output,True)

    for p in range(num_processes):
        # process communication files selection (one per process) 
        p_comm_filename = process_communication_filename + "_" + str(p)
        p_comm_filename_output = test_name + "_pruned_" + str(num_processes) + "_" + str(p)
        if os.path.isfile(p_comm_filename):
            select_rows_from_csv_irregular_columns(selected_rows,p_comm_filename,p_comm_filename_output,True)
        
        # trace files selection (one per process)

        t_filename = trace_filename + "_" + str(p) + "_trace"
        t_filename_output = test_name + "_pruned_" + str(num_processes) + "_" + str(p) + "_trace"
        if os.path.isfile(t_filename):
            select_rows_from_csv_irregular_columns(selected_rows,t_filename,t_filename_output,False)
        




