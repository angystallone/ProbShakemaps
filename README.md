# ProbShakemaps

**ProbShakemaps** is an open source Python toolbox that generates Probabilistic Shakemaps, an evolved version of traditional [Shakemaps](https://github.com/DOI-USGS/ghsc-esi-shakemap). 

Dependencies
-----------------------------

 * All [Shakemap](https://github.com/DOI-USGS/ghsc-esi-shakemap) and [OpenQuake](https://github.com/gem/oq-engine/blob/master/README.md) dependencies
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
            

REQUIRED TO RUN
------------------

1) <ins>Shakemap Docker Image</ins> --> Install it from [USGS_Shakemap_Image](https://hub.docker.com/r/seddev/shakemap). Results shown here and in the article have been generated using the [INGV Shakemap Image](https://github.com/INGV/shakemap-input-eu). The folder `data` must contain a subfolder named as the ID of the event, containing the following files: `shake_result.hdf`, `event.xml`, `stationlist.json` (optional, needed for tools 'StationRecords' and 'GetStatistics'). The Vs30 grd file must be but in the `data/shakemap_data/vs30`.
2) <ins>Ensemble file with scenarios list</ins> --> The following scenario parameterization is required: magnitude, longitude, latitude, hypocenter depth (km), strike, dip, rake, fault area (L x W, *km^2*), fault length (L, *km*), slip (*m*). The file must be put in the folder INPUT_FILES/ENSEMBLE. The last 2 characters of the ensemble filename should indicate the ensemble sample number (i.e. `_01.txt` for ensemble 1, `_02.txt` for ensemble 2...). See an example of ensemble file at [list_scenarios.txt](https://github.com/angystallone/ProbShakemaps/blob/main/INPUT_FILES/ENSEMBLE/list_scenarios_01.txt).
3) <ins>POIs file</ins> --> two space-separated columns .txt file with LAT and LON of the POIs. See an example of POIs file at [POIs_grid.txt](https://github.com/angystallone/ProbShakemaps/blob/main/INPUT_FILES/POIs_grid.txt).
4) <ins>input_file.txt</ins> --> File containing the parameters required for running the probabilistic analysis (do not rename it). See an example of input file at [input_file.txt](https://github.com/angystallone/ProbShakemaps/blob/main/INPUT_FILES/input_file.txt).
  
> * TectonicRegionType: as defined in `OpenQuake` tectonic regionalisation.
> * Magnitude_Scaling_Relationship: as required from `openquake.hazardlib.scalerel`.
> * ID_Event: `Shakemap` ID of the event.
> * POIsfile: two space-separated columns txt file with LAT and LON of the POIs.
> * Vs30file: GMT grd Vs30 file; if not provided, default Vs30 value (760 m/s) is used.
> * CorrelationModel: as required from `openquake.hazardlib.correlation`.
> * CrosscorrModel: as required from `openquake.hazardlib.cross_orrelation`.
> * vs30_clustering: `True` value means that Vs30 values are expected to show clustering (as required from `openquake.hazardlib.correlation`).
> * truncation_level: number of standard deviations for truncation of the cross-correlation model distribution (as required from `openquake.hazardlib.cross_correlation`).
> * NumGMPEsRealizations: total number of GMPEs realizations at a specific POI for spatial and crosscorrelation of ground motion fields.
> * num_processes: number of CPU cores for code parallelization.


Usage
------------------

**RUN PROBABILISTIC ANALYSIS**

```bash
python ProbShakemaps.py --task RunProbAnalysis
```

OUTPUT

`SIZE_{num_scenarios}_ENSEMBLE_{num_ensemble}_{imt}.hdf5`: HDF5 file generated by the probabilistic analysis. Each POI is assigned an IMT distribution including as many values as the number of GMPEs realizations chosen by the user.

`params.json`: File with parameters required by the ProbShakemaps tools

***************************************

**GENERATE PROBSHAKEMAPS**

Include two utility tools ('StationRecords' and 'QueryHDF5') and the main tools to generate Probabilistic Shakemaps: 'GetStatistics' and 'EnsemblePlot'. Probabilistic Shakemaps helps visualizaing ensemble predictions at the POIs.

**TOOL: 'StationRecords'**

Inspect Shakemap .json station file

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool StationRecords --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json
```
OUTPUT

`Data_stationfile_{imt}.png`: Plot data from Shakemap .json station file for the selected IMT (PGA in the example).

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/Data_stationfile_PGA.png" alt="Data_stationfile_PGA.png" width="60%" height="60%">
</p>

**TOOL: 'QueryHDF5'**

Navigate and query the .HDF5 file 

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool QueryHDF5 --imt PGA --scenario 50 --pois_file POIs.txt
```

OUTPUT

`GMF_info.txt`: Print the ground motion fields for the selected scenario at the POIs listed in `POIs.txt`.

Preview of the output file:
<pre>
GMF realizations at Site_LAT:43.2846_LON:12.7778 for Scenario_10: [0.17520797, 0.21844997, 0.093965515, 0.27266037, 0.079073295, 0.09725358, 0.08347481, 0.06693749, 0.005907976, 0.060873847]
GMF realizations at Site_LAT:43.1846_LON:12.8778 for Scenario_10: [0.100996606, 0.35003924, 0.24363522, 0.19941418, 0.15757227, 0.1009447, 0.19146584, 0.06460667, 0.03146108, 0.097111605]
GMF realizations at Site_LAT:43.0846_LON:13.4778 for Scenario_10: [0.18333985, 0.11954803, 0.2914887, 0.050770156, 0.07628956, 0.17871241, 0.10297835, 0.15162756, 0.020328628, 0.04087482]
</pre>

**TOOL: 'GetStatistics'**

* Calculates and save statistics ('Mean','Median','Percentile 10','Percentile 20','Percentile 80','Percentile 90'). If the user does not provide a file with scenarios weights, they are considered equally probable.
* Plots calculated statistics at each POI.
* Plots the IMT cumulative distribution and main statitics at a specific POI together with the estimated IMT value at the closest station (datum taken from the Shakemap .json station file). Note: the IMT cumulative distribution is based an all scenarios in the ensemble.

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool GetStatistics --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json --pois_file POIs.txt
```

OUTPUT

* npy files with statistics saved in the `npyFiles` folder
* Statitics map distributions saved in the `STATISTICS` folder

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/STATISTICS/summary_stats_forReadMe.png" alt="SummaryStats" width="85%" height="85%">
</p>

* Plot of Datum-Ensemble comparison at each POI saved in the `DISTRIBUTIONS` folder

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/DISTRIBUTIONS/Distr_POI-003.png" alt="DatumEnsemble" width="80%" height="80%">
</p>

**TOOL: 'EnsemblePlot'**

Summarize the IMT distribution at the selected POIs with a boxplot. Note: the IMT distribution is based an all scenarios in the ensemble.

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool EnsemblePlot --imt PGA --pois_file POIs.txt
```
OUTPUT

`POIs_Map.pdf`: Spatial map of the POIs

`Ensemble_Spread_Plot_{imt}.pdf`: Boxplot

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/summary_stats_forReadMe.png" alt="DatumEnsemble" width="80%" height="80%">
</p>

**POIs SUBSET OPTION**

When using the tools 'QueryHDF5', 'GetStatistics' and 'EnsemblePlot', the user can also require to extract a subset of POIs within a maximum distance from the event epicenter following one of these two possible spatial distributions: <ins>random</ins> and <ins>azimuthally uniform</ins>. This changes the command line to ('GetStatistics' example):

```bash
python ProbShakemaps.py --task GenerateProbShakemap --tool GetStatistics --imt PGA --imt_min 0.01 --imt_max 10 --station_file stationlist.json --pois_file POIs.txt --pois_subset n_pois 10 max_distance 50 --pois_selection_method azimuth_uniform
```

**HPC**

The code can be run on a cluster enjoying parallelization. See an example of bash file to run the code on a HPC cluster at [run_code.bash](https://github.com/angystallone/ProbShakemaps/blob/main/run_code.bash). IMPORTANT: the number set at `--ntasks-per-node` must coincide with `num_processes` in `input_file.txt`.


License
-----------------------------

This project is released under the [MIT License](LICENSE).


Contact
-----------------------------

If you need support, write to [angela.stallone@ingv.it](mailto:angela.stallone@ingv.it).
