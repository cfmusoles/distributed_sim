#!/bin/bash
#$ -pe openmpi-ib 4
#$ -l h_rt=08:00:00
#$ -l rmem=3G
#$ -M c.f.musoles@sheffield.ac.uk
#$ -m a

APP_PATH="./"
APP_NAME="distSim"

PROCESSES=4    #should be the same as mpi processes requested
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
module load mpi/gcc/openmpi
module load compilers/gcc

NODE_INFO="node_info_"$PROCESSES
cat $PE_HOSTFILE > $NODE_INFO

ITERATIONS=1

for i in $(seq 1 $ITERATIONS)
do
  	SEED=$RANDOM
        echo "Iteration $i"
        mpirun -np $PROCESSES $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_random" $COMM_PATTERN "random" $SEED
        sleep 20
        #mpirun -np $PROCESSES $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_randomBalanced" $COMM_PATTERN "randomBalanced" $SEED
        #sleep 20
        #mpirun -np $PROCESSES $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_graphPartitioning" $COMM_PATTERN "graphPartitioning" $SEED
        #sleep 20
        #mpirun -np $PROCESSES $TRACING $APP_PATH$APP_NAME $EXPERIMENT_NAME"_graphPartitioning_repartition" $COMM_PATTERN "graphPartitioning" $SEED $EXPERIMENT_NAME"_graphPartitioning__"$PR$
        #sleep 20
done

# mpi run debugging command: --mca oob_tcp_debug 99 --mca oob_tcp_verbose 99 --mca mca_verbose 99 --mca pml_base_verbose 99 --mca btl_base_debug 99 --mca btl_base_verbose 99 --mca mtl_base_verbose 99

# mpi run mlt: --mca mtl psm2
