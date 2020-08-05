#!/bin/bash


export CH="cd /home/bsmithers/software_dev/Analysis/plotting"
export ENV="py3v401"
export RU="python compare_plots.py -i /data/user/bsmithers/processed_runs/ -c viridis --sep_matter"

# Connect to the server, go to the plotting script, engage the environment, and run the script
ssh bsmithers@cobalt.icecube.wisc.edu "${CH} ; ${ENV} ; ${RU} ; logout"

# Copy the produced files to the current cirectory 
export WHERE_TO=${PWD}
scp bsmithers@cobalt.icecube.wisc.edu:/data/user/bsmithers/processed_runs/output/*.png ${WHERE_TO}
