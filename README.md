# ProbShakemaps

The **ProbShakemaps** is an open source Python toolbox that generates a probabilistic version of the traditional [Shakemaps](https://github.com/DOI-USGS/ghsc-esi-shakemap). 

Dependencies
-----------------------------

 * All [Shakemap](https://github.com/DOI-USGS/ghsc-esi-shakemap) dependencies
 * basemap
 * seaborn
 
 For more details, check `requirements.txt` for a list of all the installed Python packages and their versions.
 
Command line usage
------------------
<pre>
usage: ProbShakemaps.py [-h] --task {RunProbAnalysis,GenerateProbShakemap} [--imt {PGA,PGV}] 
                        [--tool {StationRecords,QueryHDF5,GetStatistics,EnsemblePlot}] [--imt_min IMT_MIN] 
                        [--imt_max IMT_MAX] [--station_file STATION_FILE] [--scenario SCENARIO] [--pois_file POIS_FILE] 
                        [--ensemble_number ENSEMBLE_NUMBER] [--deg-round DEG_ROUND] [--pois_subset] [--n_pois N_POIS]
                        [--max_distance MAX_DISTANCE] [--pois_selection_method {random,azimuth_uniform}] 
                        [--fileScenariosWeights FILESCENARIOSWEIGHTS]

ProbShakemap Toolbox

optional arguments:
  -h, --help            show this help message and exit

input params:
  --task {RunProbAnalysis,GenerateProbShakemap}
                        Task to perform
  --imt {PGA,PGV}       Intensity measure type (IMT)
  --tool {StationRecords,QueryHDF5,GetStatistics,EnsemblePlot}
                        Tool to use
  --imt_min IMT_MIN     Minimum value for the selected IMT (for plot only)
  --imt_max IMT_MAX     Maximum value for the selected IMT (for plot only)
  --station_file STATION_FILE
                        Shakemap .json station file
  --scenario SCENARIO   Scenario number
  --pois_file POIS_FILE
                        Filename with latitude and longitude of POIs
  --ensemble_number ENSEMBLE_NUMBER
                        Ensemble sample number
  --deg-round DEG_ROUND
                        Rounding precision for latitude and longitude
  --pois_subset         Extract a subset of POIs
  --n_pois N_POIS       Number of POIs in the subset
  --max_distance MAX_DISTANCE
                        Max distance from epicenter of POIs in the subset
  --pois_selection_method {random,azimuth_uniform}
                        Selection method for the POIs of the subset
  --fileScenariosWeights FILESCENARIOSWEIGHTS
                        File with scenarios weights
</pre>                        
            
            
Examples
------------------

**RUN PROBABILISTIC ANALYSIS**

```bash
python ProbShakemaps.py --task RunProbAnalysis
```

**GENERATE PROBSHAKEMAPS**

* TOOL: 'StationRecords'

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool StationRecords --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json
```

* TOOL: 'QueryHDF5'

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool QueryHDF5 --imt PGA --scenario 50 --pois_file POIs.txt
```

* TOOL: 'GetStatistics'

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool GetStatistics --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json --pois_file POIs.txt
```

* TOOL: 'EnsemblePlot'

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool EnsemblePlot --imt PGA --pois_file POIs.txt
```

License
-----------------------------

This project is released under the [MIT License](LICENSE).
