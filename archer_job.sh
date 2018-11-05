#!/bin/bash --login

# name of the job
#PBS -N simulation_job  
# how many nodes
#PBS -l select=8
# walltime
#PBS -l walltime=0:20:0
# budget code
#PBS -A e582  

APP_NAME="distSim"
ITERATIONS=3
REPETITIONS=2
PROCESSES=192
EXPERIMENT_NAME="test"

#two args:
# $1 > distribution algorithm used by simulator
# $2 > if results should be pruned
run_experiment() {
	P="$1"
	DISTRIBUTION="$2"
	COMM_PATTERN="$3"
	SEED="$4"
	for i in $(seq 1 $REPETITIONS)
	do
		aprun -n $P $APP_NAME -n $EXPERIMENT_NAME -c $COMM_PATTERN -p $DISTRIBUTION -s $SEED -k 1000 -f 200 -t 350 -m "mvc"
		sleep 1
	done
	
}




# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

for i in $(seq 1 $ITERATIONS)
do
	SEED=$RANDOM
	run_experiment $PROCESSES "roundrobin" "pex" $SEED
	run_experiment $PROCESSES "hypergraphPartitioning" "nbx" $SEED
done


#cat output*.out > helloworld.out

