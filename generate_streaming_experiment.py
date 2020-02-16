# Create ARCHER job files based on parameters passed
# snn_praw: This experiment focuses on comparing runtime impact on SNN simulations when using alternative partitioning methods to distribute workloads:
#	- zoltan
#	- hyperPRAW (architecture aware streaming)
# Both approaches are tested with and without previous neuron activity information


import sys
import math

#job templates
template_1 = '''#!/bin/bash --login

# name of the job
#PBS -N '''
template_2 = '''
# how many nodes
#PBS -l select='''
template_3=''':bigmem='''
template_4='''
# walltime
#PBS -l walltime=24:00:0
# budget code
#PBS -A e582

APP_NAME="distSim"
PRUNE_COLUMN=4 # propagation time (sync time + idle)
ITERATIONS=2
REPETITIONS=1
PROCESSES='''
template_5='''
EXPERIMENT_NAME='''
template_6='''
# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

# bandwidth probing parameters
SIZE=512
ITERATIONS=20
WINDOW=10
#renaming is necessary to avoid clashes between simultaneous jobs
ORIGINAL_BM_FILE="results_mpi_send_bandwidth_"$PROCESSES
aprun -n $PROCESSES mpi_perf $SIZE $ITERATIONS $WINDOW
for p in $(seq 1 10)
do
	FILENAME="results_mpi_send_bandwidth_"$p"_"$PROCESSES
	if [ ! -f $FILENAME ]; then
	    BM_FILE="results_mpi_send_bandwidth_"$p"_"$PROCESSES
	    break
	fi
done

mv $ORIGINAL_BM_FILE $BM_FILE


# run experiments

COMM_PATTERN="nbx"
SEED=$RANDOM

for i in $(seq 1 $REPETITIONS)
do
	# run baseline for all partitioning candidates
	aprun -n $PROCESSES $APP_NAME -n $EXPERIMENT_NAME -c $COMM_PATTERN -p "prawV_par:hypergraphPartitioning:prawV_seq" -s $SEED -k 1000 -f 160 -t 700 -m "mvc" -i 24 -b $BM_FILE -q $ITERATIONS
	sleep 1
	
	# run with neuron activity info (only supported by prawV)
	aprun -n $PROCESSES $APP_NAME -n $EXPERIMENT_NAME"_neuron_activity" -c $COMM_PATTERN -p "prawV_par:prawV_seq" -s $SEED -k 1000 -f 160 -t 700 -m "mvc" -i 24 -b $BM_FILE -W -q $ITERATIONS
	sleep 1

done


'''


if len(sys.argv) < 6:
	print("Input error: usage -> python generate_archer_job.py jobName min_processes num_experiments geometric_step big_mem[true|false]")
	exit()

test_name = sys.argv[1]
min_processes = int(sys.argv[2])
num_experiments = int(sys.argv[3])
geometric_step = int(sys.argv[4])
big_mem = (sys.argv[5] == "true" or sys.argv[5] == "True")

process_counts = [min_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
print("Generating experiments")
print(process_counts)

for p in process_counts:
	nodes = max(int(math.ceil(p / 24)),1)
	writebuffer = open("archer_job_" + test_name + "_" + str(p) + ".sh",'w')
	writebuffer.write(template_1 + test_name)
	writebuffer.write(template_2 + str(nodes))
	writebuffer.write(template_3 + str(big_mem).lower())
	writebuffer.write(template_4 + str(p))
	writebuffer.write(template_5 + test_name)
	writebuffer.write(template_6)



