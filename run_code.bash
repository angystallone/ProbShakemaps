#!/bin/bash
#SBATCH --job-name=ProbShakemaps
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=job.out
#SBATCH --error=job.err

# Set up the environment
export OMP_NUM_THREADS=1

#######################################################

# Run Probabilistic Analysis
python ProbShakemaps.py --task RunProbAnalysis

#######################################################

# Generate ProbShakemaps

# TOOL: 'StationRecords'
python ProbShakemaps.py --task GenerateProbShakemap --tool StationRecords --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json

# TOOL: 'QueryHDF5'
python ProbShakemaps.py --task GenerateProbShakemap --tool QueryHDF5 --imt PGA --scenario 50 --pois_file POIs.txt

# TOOL: 'GetStatistics'
python ProbShakemaps.py --task GenerateProbShakemap --tool GetStatistics --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json --pois_file POIs.txt

# TOOL: 'EnsemblePlot'
python ProbShakemaps.py --task GenerateProbShakemap --tool EnsemblePlot --imt PGA --pois_file POIs.txt