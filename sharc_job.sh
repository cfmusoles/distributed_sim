#!/bin/bash
#$ -pe mpi 48
#$ -l h_rt=08:00:00
#$ -l rmem=4G
#$ -l excl=false
#$ -M c.f.musoles@sheffield.ac.uk
#$ -m a

APP_PATH="./"
APP_NAME="distSim"

PROCESSES=48    #should be the same as mpi processes requested
EXPERIMENT_NAME="$1"
COMM_PATTERN="$2"

if [ "$3" != "" ]
then
	TRACE_NAME="$3"
	TRACING=$APP_PATH$TRACE_NAME
else
	TRACING=""
fi

#load modules
module load mpi/openmpi/2.1.1/gcc-5.4
module load dev/gcc/5.4

NODE_INFO="node_info_"$PROCESSES
cat $PE_HOSTFILE > $NODE_INFO

ITERATIONS=5

for i in $(seq 1 $ITERATIONS)
do
	SEED=$RANDOM
        echo "Iteration $i"
        mpirun -np $PROCESSES --map-by ppr:1:core $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_random" $COMM_PATTERN "random" $SEED
        sleep 20
        mpirun -np $PROCESSES --map-by ppr:1:core $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_randomBalanced" $COMM_PATTERN "randomBalanced" $SEED
        sleep 20
        mpirun -np $PROCESSES --map-by ppr:1:core $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_graphPartitioning" $COMM_PATTERN "graphPartitioning" $SEED
        sleep 20
	mpirun -np $PROCESSES --map-by ppr:1:core $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_graphPartitioning_repartition" $COMM_PATTERN "graphPartitioning" $SEED $EXPERIMENT_NAME"_graphPartitioning__"$PROCESSES"_neuron_activity"
        sleep 20
done

# mpi run debugging command: --mca oob_tcp_debug 99 --mca oob_tcp_verbose 99 --mca mca_verbose 99 --mca pml_base_verbose 99 --mca btl_base_debug 99 --mca btl_base_verbose 99 --mca mtl_base_verbose 99
