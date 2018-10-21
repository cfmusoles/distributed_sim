APP_NAME="distSim"

ITERATIONS=15

MIN_PROCESSES="$1"
MAX_PROCESSES="$2"
MODEL="$3"

run_experiments() {
	COMM_PATTERN="$1"
	PARTITIONING="$2"
	RESULTS=$MODEL"_"$COMM_PATTERN"_"$PARTITIONING
	for i in $(seq 1 $ITERATIONS)
	do
		for j in $(seq 2 $PROCESSES)
		do
			echo "Iteration $i"
			mpirun -np $j $APP_NAME -n $RESULTS -c $COMM_PATTERN -p $PARTITIONING
			sleep 1
		done
	done
}

run_partitioning_experiment() {
	MODEL_NAME="$1"
	COMM_PATTERN="$2"
	PARTITIONING="$3"
	SEED="$4"
	if [ $5 == "yes" ]
	then
		RESULTS=$MODEL_NAME"_"$COMM_PATTERN"_"$PARTITIONING"_repartition"
	else
		RESULTS=$MODEL_NAME"_"$COMM_PATTERN"_"$PARTITIONING
	fi
	
	for j in $(seq $MIN_PROCESSES $MAX_PROCESSES)
	do
		echo "Iteration $i"
		if [ $5 == "yes" ]
		then
			ACTIVITY_FILE=$MODEL_NAME"_"$COMM_PATTERN"_"$PARTITIONING"__"$j"_neuron_activity"
		else
			ACTIVITY_FILE=""
		fi
		mpirun -np $j $APP_NAME -n $RESULTS -c $COMM_PATTERN -p $PARTITIONING -s $SEED -w $ACTIVITY_FILE
		sleep 1
	done
}


## RANDOM PARTITIONING
#run_experiments "allGather" "random"
#run_experiments "p2p" "random"
#run_experiments "subscriber" "random"

## GRAPH PARTITIONING
#run_experiments "allGather" "graphPartitioning"
#run_experiments "p2p" "graphPartitioning"
#run_experiments "subscriber" "graphPartitioning"


#TO TEST REPARTITIONING

 
for i in $(seq 1 $ITERATIONS)
do
	RND=$RANDOM
	run_partitioning_experiment $MODEL "subscriber" "random" $RND "no"
	#run_partitioning_experiment $MODEL "subscriber" "graphPartitioning" $RND "no"
	#run_partitioning_experiment $MODEL "subscriber" "graphPartitioning" $RND "yes"
done
