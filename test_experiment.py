# Create ARCHER job files based on parameters passed

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
#PBS -l walltime=6:00:0
# budget code
#PBS -A e582

APP_NAME="distSim"
TEST_ITERATIONS=1
TEST_REPETITIONS=2
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
# bandwidth matrix creation
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


#two args:
# $1 > distribution algorithm used by simulator
# $2 > if results should be pruned
run_experiment() {
	P="$1"
	DISTRIBUTION="$2"
	COMM_PATTERN="$3"
	SEED="$4"
	for i in $(seq 1 $TEST_REPETITIONS)
	do
		#aprun -n $P $APP_NAME -n $EXPERIMENT_NAME -c $COMM_PATTERN -p $DISTRIBUTION -s $SEED -k 500 -f 1000 -t 150 -m "cm" -i 24 -b $BM_FILE
		sleep 1
		aprun -n $P $APP_NAME -n $EXPERIMENT_NAME -c $COMM_PATTERN -p $DISTRIBUTION -s $SEED -k 80 -f 100 -t 300 -m "mvc" -i 24 -b $BM_FILE
		sleep 1
	done
}

for r in $(seq 1 $TEST_ITERATIONS)
do
	SEED=$RANDOM
	#run_experiment $PROCESSES "roundrobin" "nbx" $SEED 
	run_experiment $PROCESSES "prawE" "nbx" $SEED 
	run_experiment $PROCESSES "prawE_without" "nbx" $SEED 
	run_experiment $PROCESSES "hypergraphPartitioning" "nbx" $SEED 
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



