
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
ITERATIONS=2
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
        HGRAPH="$5"
	if [ "$6" != "" ]
        then
                BANDWIDTH_FILE=""
                FILENAME=$EXPERIMENT_NAME"_"$DISTRIBUTION"_noB"
        else
                BANDWIDTH_FILE="-b ""$6"       
                FILENAME=$EXPERIMENT_NAME"_"$DISTRIBUTION"_yesB"

        fi

        aprun -n $P $APP_NAME -n $FILENAME -c $COMM_PATTERN -p $DISTRIBUTION -s $SEED -h $HGRAPH -N $BANDWIDTH_FILE
        sleep 1
}




# This shifts to the directory that you submitted the job from
cd $PBS_O_WORKDIR

# generate comm bandwidth files for current architecture
aprun -n $PROCESSES mpi_perf

for i in $(seq 1 $ITERATIONS)
do
        SEED=$RANDOM
        run_experiment $PROCESSES "random" "nbx" $SEED "192bit.2.mtx.hgr" "results_mpi_send_bandwidth_"$PROCESSES
        run_experiment $PROCESSES "zoltanFile" "nbx" $SEED "192bit.2.mtx.hgr" "results_mpi_send_bandwidth_"$PROCESSES
        run_experiment $PROCESSES "praw" "nbx" $SEED "192bit.2.mtx.hgr" "results_mpi_send_bandwidth_"$PROCESSES
done

