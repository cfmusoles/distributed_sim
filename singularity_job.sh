#!/bin/bash
PROCESSES=48    #should be the same as mpi processes requested
#$ -pe mpi 48
#$ -l h_rt=01:00:00
#$ -l rmem=2G
# Email notifications
#$ -M c.f.musoles@sheffield.ac.uk
# Email notifications if the job aborts
#$ -m a

IMAGE_PATH="/home/cop12c/"
IMAGE_NAME="distributed_sim.img"
APP_PATH="/home/cop12c/distributed_sim/"
APP_NAME="distSim"

#load modules
module load mpi/openmpi/2.0.1/gcc-4.9.4

NODE_INFO = "node_info_"$PROCESSES
cat $PE_HOSTFILE > $NODE_INFO

ITERATIONS=10

# compile code using singularity image
#cd $APP_PATH
#singularity exec $IMAGE_PATH$IMAGE_NAME make
#cd ~
#chmod +x $APP_PATH$APP_NAME


for i in $(seq 1 $ITERATIONS)
do
        echo "Iteration $i"
        mpirun -np $PROCESSES singularity exec $IMAGE_PATH$IMAGE_NAME $APP_PATH$APP_NAME results
        sleep 20
done
