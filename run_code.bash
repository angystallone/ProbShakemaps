#!/bin/bash
#SBATCH --job-name=ProbShakemap
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=job.out
#SBATCH --error=job.err

# Set up the environment
export OMP_NUM_THREADS=1

#######################################################

# Tools 

# TOOL: 'StationRecords'
python ProbShakemap.py --imt PGA --tool StationRecords --imt_min 0.01 --imt_max 10 --station_file stationlist.json

# TOOL: 'Save_Output'
python ProbShakemap.py --imt PGA --tool Save_Output --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10

# TOOL: 'QueryHDF5'
python ProbShakemap.py --imt PGA --tool QueryHDF5 --scenario 50 --pois_file POIs.txt


#######################################################

# Prob_tools

# TOOL: 'GetStatistics'
python ProbShakemap.py --imt PGA --prob_tool GetStatistics --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 1

# TOOL: 'GetDistributions'
python ProbShakemap.py --imt PGA --prob_tool GetDistributions --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 10 --station_file stationlist.json

# TOOL: 'EnsemblePlot'
python ProbShakemap.py --imt PGA --prob_tool EnsemblePlot --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10
