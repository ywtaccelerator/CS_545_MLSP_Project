#!/bin/bash
#SBATCH --account=def-mfrase9
#SBATCH --time=0-4:30:00 # maximum running time [format: days-hours:minutes:seconds]
#SBATCH --gres=gpu:v100l:1 # number of GPUs will be used
#SBATCH --cpus-per-task=4 # number of CPUs will be used
#SBATCH --mem=95990M # size of memory will be used
#SBATCH --nodes=1
#SBATCH --job-name=run_example_python_test_01
#SBATCH --output=log_%x_jobID-%j.out


module load nixpkgs/16.09
module load python/3.7.4
module load scipy-stack/2019b
module load gcc/7.3.0
module load cuda/10.1
module load cudacore/.10.1.243
module load cudnn/7.6.5

source ~/projects/def-mfrase9/wyang010/for_tmpdir_use/ak5/bin/activate
pip install --no-index --upgrade pip

python *.py
