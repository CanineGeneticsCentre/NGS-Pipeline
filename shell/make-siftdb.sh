#!/usr/bin/env bash

# We assume running this from the script directory
job_directory="${HOME}/hpc-work/SIFT/"
#data_dir="${HOME}/hpc-work/SIFT/sift-test_files/"

job_file="${job_directory}/make-siftdb.job"

echo "#!/usr/bin/env bash
#SBATCH -A MELLERSH-SL3-CPU
#SBATCH --output=/home/%u/hpc-work/logs/job-%j.out
#SBATCH --error=/home/%u/hpc-work/logs/job-%j..err
#SBATCH -p icelake
#SBATCH --time=00:01:00
#SBATCH --job-name=sift.job
#SBATCH --mail-type=ALL

source $CONDA_PREFIX/etc/profile.d/conda.sh  # Always add this command to your scripts
conda activate ENA

perl /rds/project/rds-Qr3fy2NTCy0/Software/local/src/make-SIFT-db-all.pl -config ${HOME}/hpc-work/SIFT/sift-test_files/homo_sapiens-test.txt" > $job_file

#sbatch $job_file
echo sbatch $job_file;