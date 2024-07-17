#!/bin/bash
#SBATCH --job-name=probshakemap
#SBATCH -p big
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=00:30:00
#SBATCH --output=job.out
#SBATCH --error=job.err

export OMP_NUM_THREADS=1
export PATH=/shakemap/home/shake/miniconda/envs/shakemap/bin/:/home/shake/miniconda/envs/shakemap/bin/:/home/shake/miniconda/envs/env/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

module load singularity/3.4.1

singularity exec -B $(pwd)/..:/home/shake/shakemap_profiles -B $(pwd)/../../shakemap_data:/home/shake/shakemap_data -B $(pwd)/../../../data/local:/home/shake/.local ../../../../shakemap/myimage_new.sif /home/shake/miniconda/envs/shakemap/bin/python3 /home/shake/shakemap_profiles/world/ProbShakemap.py --imt PGA --prob_tool GetStatistics --num_processes 8 --pois_file grid_10km_no_stations.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 1 --deg_round 2 --vector_npy
# singularity exec -B $(pwd)/..:/home/shake/shakemap_profiles -B $(pwd)/../../shakemap_data:/home/shake/shakemap_data -B $(pwd)/../../../data/local:/home/shake/.local ../../../../shakemap/myimage_new.sif /home/shake/miniconda/envs/shakemap/bin/python3 /home/shake/shakemap_profiles/world/ProbShakemap.py --imt PGA --prob_tool GetDistributions EnsemblePlot --num_processes 8 --pois_file grid_10km_new.txt --numGMPEsRealizations 10 --reuse_pois_subset --imt_min 0.001 --imt_max 1 --deg_round 5 --station_file stationlist.json 

module unload singularity/3.4.1
