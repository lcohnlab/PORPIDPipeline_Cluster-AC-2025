# Notes
This version of the pipeline switched to using a maxfs filter to exclude artifactual families of smaller sizes. This works by selecting the largest likely_real family in the sample (before postproc filtering) and then setting a lower family size limit based on that. So a 0.1 maxfs filter would set a min fs of 10 if largest family was 100 and any families below 10 would be discarded. 

This version of the pipeline was used to run MCA0835 batch 2 and three samples from SCOPE_gates_ATI
   SCOPE1498-2022.09.21_PLA_pblib
   SCOPE2153-2022.09.23_PLA_pblib
   SCOPE2929-2022.01.31_PLA_pblib

Because there were no differences in the core PORPID code between the 50x filter branch used to run the previous MCA samples, and this branch, and all artifact filtering was done afterwards using custom R scripts, this branch should be used for publication of the MCA datasets. 

# PORPIDpipeline_FH_Cluster
FredHutch Rhino Cluster Install of PORPIDpipeline

This install and scripts are edited to run on the FH cluster set up with Rhino. The original code exists in the MurrellGroup GitHub Page: **https://github.com/MurrellGroup/PORPIDpipeline**


# How To Run PORPID Pipeline

This page contains information regarding how to run the PORPIDpipeline. There are 3 main contents:

1. Running Helper Files
2. Running PORPID via Snakemake
3. Creating Panel Files (Not included here)
Helper files help set up a run, and Panel files are necessary for PostProc in a PORPID run. The contents of this page should help you set up, and run the pipeline. Below is the directory set up of PORPID:

```
/fh/fast/cohn_l/grp/Pipelines
├── setup_data.sh
├── PORPIDpipeline
│   └── ...
└── data
   ├── in
   │   ├── config
   │   │   ├── *.yaml
   │   │   └── templates
   │   │       ├── *.csv
   │   │       └── primersheet2config_panel.py
   │   └── *.fastq.gz
   └── out
```

# Running Helper Files
To help set up PORPID runs and clean up past runs, I have made some helper files located in the Pipelines directory that can be run with the appropriate inputs from command line via FH cluster. Below are instructions to run these helper files:

## primersheet2config_panel.py

Located in the directory './PORPIDpipeline/data/in/config/templates'. Run this python file to generate $DS-config.yaml file from a given csv. Should be run from the templates directory, and needs to have python module loaded.

Run: `python primersheet2config_panel.py [filename]` (Example:  python primersheet2config_panel.py demo.csv)

## setup_data.sh

This script helps set up a new PORPID run by copying the config and raw read data from `./Pipelines/data/in/` into `./Pipelines/PORPIDpipeline`. If `config.yaml` exists already, then it will be overwritten (but also `$DS-config.yaml` copy will be kept, where $DS is the dataset name). Once set, a message will be printed to update Panel files if necessary, and a command to run the snakemake.

Run: `./setup_dataset.sh [dataset-filename]` (Example: ./set_dataset.sh demo)

# Running PORPID via Snakemake

We run PORPID via snakemake through the command `sbatch`. There are more specific ways to batch a command, but the default sbatch should work fine. The file to run the snakefile is a shell script `runPipelineCluster.sh`. 

## snakefile

The snakefile contains the rules for PORPID and also parameters that may need to be changed depending on the run. For example, with some previous runs, the max_length parameter was changed to 9000 (from the default 4300). This file mostly remains the same, other than any specific snakefile parameter changes. 

## config.yaml

This config.yaml file is the main config file used per PORPID run, dictating what samples to read from the raw-read and also the cDNA/sec_str primers and panel files. This is created from the primersheet2config.py file, and moved to the PORPID folder using the setup.sh helper files. Double check this config file looks correct before a run to ensure primers are correct/comment out any samples you want to run seperately

## Modules for Rhino

To run the PORPID pipeline on the FH cluster, you need to load the `snakemake` module. Once you SSH into rhino (ssh user-name@rhino), you can use the command `ml snakemake` to load the snakemake module. You can verify this module is loaded by using the command `module list`. (Note: this module also loads Python 3.9, which is used for the helper file primersheet2config.py)

I also suggest using TMUX to run PORPID. TMUX is a terminal multiplexer, and can ensure PORPID runs without interrupting standard terminal usage. If needed you can also grab a node (grabnode) for additional computing power, but that is not always necessary). Load tmux at the start of your rhino session by using `ml tmux`. You can view tmux sessions with the command `tmux ls` and detach from a session using `ctrl+b d`. Detaching keeps the session running, but if you want to close an active tmux session, type `exit`. To attach to your previous tmux session you can use the command `tmux attach-session`.

### The suggested workflow to run PORPID is as follows:

1. SSH into Rhino and navigate to Pipelines directory in fh/fast/cohn_l/grp
2. Load modules for TMUX/snakemake
3. Use the primersheet2config.py to create a config file for your $DS from the data/in/config/template directory
4. Back in the Pipeline directory, use ./setup $DS to set up your run in the PORPIDpipeline folder
5. Run `sbatchrunPipelineCluster.sh $DS` to run PORPID. You can view the output as a slurm-*.out log.

Once your run is complete, you can move the output files to the data/out directory
### Example Sbatch Command:
`sbatch runPipelineCluster.sh [dataset-filename]` 

*Example: sbatch runPipelineCluster.sh 2023-09_FH05_bc1009*

# Clearing Metadata

If a new user runs the pipeline after another user previously used the pipeline, you may run into a permission error such as the following: PermissionError: [Errno 13] Permission denied: '/fh/fast/cohn_l/grp/Pipelines/sga-index-consensus/.snakemake/metadata. In such a case you would need to clear the metadata directory in the hidden snakemake directory. 

Run the command rm -rf .snakemake/metadata/* to clear out the metadata. (Note: For best practice, at the end of a complete run, you should clear out the metadata for future users). 
