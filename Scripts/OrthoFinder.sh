#!/bin/bash
#!/bin/bash
#SBATCH --job-name=ortho             # job name
#SBATCH --ntasks=10                  # number of tasks across all nodes
#SBATCH --partition=general          # name of partition to submit job
#SBATCH --time=8-00:00:00
#SBATCH --output=ortho_run3.out          # Output file. %j is replaced with job ID
#SBATCH --error=ortho_run3.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=azk0151@auburn.edu
     #
     #  load the module
     module load orthofinder/2.4.0     #
     # OMP_NUM_THREADS much match the number of cores requested from the queue
   orthofinder -f ./renamed_faa_files -t 8

