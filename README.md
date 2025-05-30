# How to run

## First create the env and activate it
'''sh
mamba env create -f envs/environment.yml
mamba activate rna_editing_env
'''

## Next, run snakefile to start the analysis
'''sh
snakemake -j 16 --printshellcmds -n # Dry run
snakemake -j 16 --printshellcmds --latency-wait 30 # You can change the number of tasks running parallel (-j)
'''