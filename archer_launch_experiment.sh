#!/bin/bash

APP_NAME="distSim"
EXPERIMENT_NAME="$1"
MIN_PROCESSORS="$2"
NUM_EXPERIMENTS="$3"
GEOMETRIC_STEP="$4"
BIG_MEM_NODE="$5"
#create working directory
mkdir $WORK_DIR/$EXPERIMENT_NAME

#generate job files
python generate_archer_job.py $EXPERIMENT_NAME $MIN_PROCESSORS $NUM_EXPERIMENTS $GEOMETRIC_STEP $BIG_MEM_NODE
sleep 1

#geometric progression
for p in $(seq 1 $NUM_EXPERIMENTS)
do
	PROCESS_COUNT=$(($MIN_PROCESSORS * $GEOMETRIC_STEP ** ($p-1)))
	FILENAME="archer_job_"$EXPERIMENT_NAME"_"$PROCESS_COUNT".sh"
	cp $FILENAME $WORK_DIR/$EXPERIMENT_NAME/
	rm $FILENAME
done

#copy necessary files
cp $APP_NAME $WORK_DIR/$EXPERIMENT_NAME/
cp resources/multiareaModel.json $WORK_DIR/$EXPERIMENT_NAME/
chmod +x archer_retrieve_results.sh
cp archer_retrieve_results.sh $WORK_DIR/$EXPERIMENT_NAME/
cp prune_results.py $WORK_DIR/$EXPERIMENT_NAME/

#submit job files
cd $WORK_DIR/$EXPERIMENT_NAME

#geometric progression
for p in $(seq 1 $NUM_EXPERIMENTS)
do
	PROCESS_COUNT=$(($MIN_PROCESSORS * $GEOMETRIC_STEP ** ($p-1)))
	FILENAME="archer_job_"$EXPERIMENT_NAME"_"$PROCESS_COUNT".sh"
	echo "Launching job: "$FILENAME
	# launching trice to get three different node allocations
	qsub $FILENAME
	qsub $FILENAME
	rm $FILENAME
done

