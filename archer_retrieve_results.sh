#!/bin/bash

EXPERIMENT_NAME="$1"
APP_DATA_FOLDER="/home/e582/e582/musoles/projects/results"

tar -zcvf $EXPERIMENT_NAME".tar.gz" $EXPERIMENT_NAME*
cp $EXPERIMENT_NAME".tar.gz" $APP_DATA_FOLDER/
cd ..
rm -r $EXPERIMENT_NAME
cd $APP_DATA_FOLDER
