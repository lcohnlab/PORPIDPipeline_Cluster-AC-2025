#!/bin/bash
# Creates a new dataset to run PORPIDpipeline on in the user's home directory
# run from /fh/fast/cohn_l/grp/Pipelines

# if no arguments were passed, display usage info and run list_datasets
if [ $# -eq 0 ]; then
  echo "Usage: $0 <dataset name>"
  echo "For a list of datasets, run: list_datasets"
  exit 1
fi

DS=$1 # assign DS variable to first arg

# Check if the current working directory matches the expected path
if [ "$(pwd)" != "/fh/fast/cohn_l/grp/Pipelines" ]; then
  echo "This script should be executed from the directory '/fh/fast/cohn_l/grp/Pipelines'"
  exit 1
fi

# make sure dataset doesn't have odd characters or spaces

# # make sure ~/PORPIDpipeline does not already exist, otherwise print warning and exit
# if [ -d /fh/fast/cohn_l/grp/Pipelines/PORPIDpipeline/ ]; then
#   echo "./PORPIDpipeline/ already exists, canceling dataset creation. Delete it with: delete_dataset $DS"
#   exit 1
# fi

# # make sure ~/DS does not already exist, otherwise print warning and exit
# if [ -d /fh/fast/cohn_l/grp/Pipelines/$DS ]; then
#   echo "$DS already exists, canceling dataset creation. Delete it with: delete_dataset $DS"
#   exit 1
# fi

# make sure /data/in/<dataset>.fastq.gz exists as does /data/in/<dataset>-config.yaml
if [ ! -f ./data/in/$DS.fastq.gz ]; then
  echo "/data/in/$DS.fastq.gz does not exist, canceling dataset creation."
  exit 1
fi

if [ ! -f ./data/in/config/$DS-config.yaml ]; then
  echo "/data/in/config/$DS-config.yaml does not exist, canceling dataset creation."
  exit 1
fi

echo "Setting up dataset $DS in your PORPIDpipeline directory..."

# # copy PORPIDpipeline to user's home directory and rename to <dataset>
# cp -r ./src/PORPIDpipeline/ ./ && mv ./PORPIDpipeline/ ./$DS #&& chgrp -R h705 ~/$DS

# copy yaml and fastq.gz files
cp ./data/in/config/$DS-config.yaml ./PORPIDpipeline/ && cp ./PORPIDpipeline/$DS-config.yaml ./PORPIDpipeline/config.yaml
cp ./data/in/$DS.fastq.gz  ./PORPIDpipeline/raw-reads/

echo "Dataset $DS has moved over to PORPIDpipeline"
echo " "
echo "Check that you have the appropriate panel files and Snakemake parameters!"
echo "Once ready, run 'sbatch runPipelineCluster.sh $DS' to begin the analysis"
