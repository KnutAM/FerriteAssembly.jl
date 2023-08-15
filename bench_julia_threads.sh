#!/bin/bash
#SBATCH --partition=standard  # Partition (queue), some alt: standard, shortrun_large, shortrun_small, testing, fat, gpu01queue, gpu02queue, gpu03queue
#SBATCH -N 1                  # Number of nodes
#SBATCH -n 20                 # Number of processes
#SBATCH --time=0-01:00:00     # Time (1h default)
# Start a job by calling 
# `sbatch <sbatch_options> <path/to/run_julia.sh> -s <path/to/julia/script.jl>`
# from the folder containing the Project.toml
# The sbatch options are any options available for the currently used sbatch, and will override settings in this file
# An example call could be 
# sbatch -J myjob -n 2 --time=2-00:00:00 $SBATCH_SCRIPTS/run_julia.sh -s runsim.jl
# Assuming that the environment variable $SBATCH_SCRIPTS contains the path to the `sbatch` folder
# Here, the job name is set to "myjob", the number of processors to 2, and the time limit to 2 days


# Get input options
while getopts ":s:" opt; do
  case $opt in
    s)
      script=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z $script ]; then
    echo "No script file given, this is mandatory!"
    exit 1
fi

# Build environment
julia --project=. --threads $SLURM_NPROCS -e 'using Pkg; Pkg.instantiate()'

# Run benchmark
echo "Threads = 1"
julia --project=. --threads 1 $script

echo "Threads = 2"
julia --project=. --threads 2 $script

echo "Threads = 4"
julia --project=. --threads 4 $script

echo "Threads = 8"
julia --project=. --threads 8 $script

echo "Threads = 16"
julia --project=. --threads 16 $script