#!/bin/bash
#SBATCH --job-name=TW3_fit            # Job name
#SBATCH --output=output.txt           # Standard output file
#SBATCH --error=error.txt             # Standard error file
#SBATCH --array=1-100
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=8             # Number of CPU cores per task
#SBATCH --time=24:00:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=alexeyvl@ucm.esm    # Email address for notifications
#SBATCH                              # This is an empty line to separate Slurm directives from the job commands

setenv OMP_NUM_THREADS 8
python TMD+tw3_replicas_cluster.py