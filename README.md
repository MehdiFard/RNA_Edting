# How to run

## First create the env and activate it
mamba env create -f envs/environment.yml
mamba activate rna_editing_env

## Next "cd" to "input_prep/" and run the snakefile to make the annotation files
snakemake -s input_prep.smk -j 16 -n # Dry run
snakemake -s input_prep.smk -j 16 --latency-wait 30 # You can change the number of tasks running parallel (-j)

## To run the analysis, "cd" to the main folder and run snakefile.
snakemake -j 16 --printshellcmds -n # Dry run
snakemake -j 16 --printshellcmds --latency-wait 30