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
ITERATIONS=1
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
	FILENAME=$EXPERIMENT_NAME"_"$DISTRIBUTION"_"$COMM_PATTERN
	aprun -n $P $APP_NAME -n $FILENAME -c $COMM_PATTERN -p $DISTRIBUTION -s $SEED
	sleep 1
}




# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

for i in $(seq 1 $ITERATIONS)
do
	SEED=$RANDOM
	run_experiment $PROCESSES "randomBalanced" "pex" $SEED
	run_experiment $PROCESSES "randomBalanced" "nbx" $SEED
done


#cat output*.out > helloworld.out

