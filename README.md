# ProbShakemaps

```ProbShakemaps``` is a Python toolbox that generates Probabilistic Shakemaps, an evolved version of traditional [Shakemaps](https://github.com/DOI-USGS/ghsc-esi-shakemap). It efficiently quantifies and propagates earthquake source uncertainty while accounting for model (GMMs) and aleatoric (GMMs variability) uncertainties. Designed for Urgent Computing applications.

Dependencies
-----------------------------

 * [Shakemap](https://github.com/DOI-USGS/ghsc-esi-shakemap) and [OpenQuake](https://github.com/gem/oq-engine/blob/master/README.md) dependencies
 * basemap
 * seaborn
 
 For more details, check `requirements.txt` for a list of all the installed Python packages and their versions.
 
Command line usage
------------------
<pre>
usage: ProbShakemaps.py [-h] [--imt {PGA,PGV}] [--tool {Save_Output,StationRecords,QueryHDF5}]
                        [--prob_tool {GetStatistics,GetDistributions,EnsemblePlot} [{GetStatistics,GetDistributions,EnsemblePlot} ...]] [--numGMPEsRealizations NUMGMPESREALIZATIONS]
                        [--num_processes NUM_PROCESSES] [--imt_min IMT_MIN] [--imt_max IMT_MAX] [--station_file STATION_FILE] [--scenario SCENARIO] [--pois_file POIS_FILE]
                        [--deg_round DEG_ROUND] [--pois_subset] [--n_pois N_POIS] [--max_distance MAX_DISTANCE] [--pois_selection_method {random,azimuth_uniform}]
                        [--fileScenariosWeights FILESCENARIOSWEIGHTS]

ProbShakemap Toolbox

optional arguments:
  -h, --help            show this help message and exit

input params:
  --imt {PGA,PGV}       Intensity measure type (IMT)
  --tool {StationRecords,Save_Output,QueryHDF5}
                        Tool(s) to use
  --prob_tool {GetStatistics,GetDistributions,EnsemblePlot} [{GetStatistics,GetDistributions,EnsemblePlot} ...]
                        ProbShakemaps Tool(s) to use
  --numGMPEsRealizations NUMGMPESREALIZATIONS
                        Total number of GMPEs random samples
  --num_processes NUM_PROCESSES
                        Number of CPU cores for code parallelization
  --imt_min IMT_MIN     Minimum value for the selected IMT (for plot only)
  --imt_max IMT_MAX     Maximum value for the selected IMT (for plot only)
  --station_file STATION_FILE
                        Shakemap .json station file
  --scenario SCENARIO   Scenario number
  --pois_file POIS_FILE
                        Filename with latitude and longitude of POIs
  --deg_round DEG_ROUND
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

1) <ins>Shakemap Docker Image</ins> --> Install it from [USGS_Shakemap_Image](https://hub.docker.com/r/seddev/shakemap). Results shown here and in the article have been generated using the [INGV Shakemap Image](https://github.com/INGV/shakemap-input-eu). The folder `data` must contain a subfolder named as the ID of the event, containing the following files: `event.xml`, `stationlist.json` (optional, needed for tools 'StationRecords' and 'GetDistributions'). The Vs30 grd file must be but in the `data/shakemap_data/vs30` (optional). NOTE: The `event.xml` file must be provided in the format required by ShakeMap (do not rename it). See an example of `event.xml` file at [event.xml](https://github.com/INGV/ProbShakemap/blob/main/event.xml). 
2) <ins>Ensemble file with scenarios list</ins> --> The following scenario parameterization is required: magnitude, longitude, latitude, hypocenter depth (km), strike, dip, rake, fault area (L x W, *km^2*), fault length (L, *km*), slip (*m*). The file must be put in the folder INPUT_FILES/ENSEMBLE. See an example of ensemble file at [list_scenarios.txt](https://github.com/angystallone/ProbShakemaps/blob/main/INPUT_FILES/ENSEMBLE/list_scenarios_01.txt).
3) <ins>POIs file</ins> --> two space-separated columns .txt file with LAT and LON of the POIs. The file must be put in the folder INPUT_FILES. See an example of POIs file at [POIs_grid.txt](https://github.com/angystallone/ProbShakemaps/blob/main/INPUT_FILES/POIs_grid.txt).
4) <ins>input_file.txt</ins> --> (Do not rename it!) File containing the parameters required for the probabilistic analysis. The file must be put in the folder INPUT_FILES. See an example of input file at [input_file.txt](https://github.com/angystallone/ProbShakemaps/blob/main/INPUT_FILES/input_file.txt).
5) <ins>fileScenariosWeights.txt</ins> --> File with scenarios weights (optional). The file must be put in the folder INPUT_FILES.

NOTE: ```ProbShakemaps``` will rely on the files provided in the folder INPUT_FILES. To run it for different events, simply rename the old INPUT_FILES folder and populate a new one from scratch.
  
> * TectonicRegionType: as defined in OpenQuake tectonic regionalisation.
> * Magnitude_Scaling_Relationship: as required from openquake.hazardlib.scalerel.
> * Rupture_aratio: rupture aspect ratio as required from openquake.hazardlib.geo.surface.PlanarSurface.from_hypocenter
> * ID_Event: Shakemap ID of the event.
> * Vs30file: GMT .grd Vs30 file; if not provided, set it to None. Default Vs30 value (760 m/s) will be used instead.
> * CorrelationModel: as required from openquake.hazardlib.correlation.
> * CrosscorrModel: as required from openquake.hazardlib.cross_orrelation.
> * vs30_clustering: `True` value means that Vs30 values are expected to show clustering (as required from openquake.hazardlib.correlation).
> * truncation_level: number of standard deviations for truncation of the cross-correlation model distribution (as required from openquake.hazardlib.cross_correlation).


Usage
------------------

**TOOLS**

`ProbShakemaps` comes with three utility tools: 'StationRecords', 'Save_Output' and 'QueryHDF5'.

**TOOL: 'StationRecords'**

Inspect Shakemap .json station file.

```bash
python ProbShakemaps.py --imt PGA --tool StationRecords --imt_min 0.01 --imt_max 10 --station_file stationlist.json
```
OUTPUT

`Data_stationfile_{imt}.pdf`: Plot data from Shakemap .json station file for the selected IMT (PGA in the example).

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/Data_stationfile_PGA.png" alt="Data_stationfile_PGA" width="60%" height="60%">
</p>


**TOOL: 'Save_Output'**

Run the probabilistic analysis and save the output to a .HDF5 file with the following hierarchical structure: 

scenario --> POI --> GMPEs realizations.

```bash
python ProbShakemaps.py --imt PGA --tool Save_Output --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10
```

OUTPUT

`SIZE_{num_scenarios}_ENSEMBLE_{IMT}.hdf5`


**TOOL: 'QueryHDF5'**

Navigate and query the .HDF5 file.

```bash
python ProbShakemaps.py --tool QueryHDF5 --imt PGA --scenario 10 --pois_file POIs.txt
```

OUTPUT

`GMF_info.txt`: Print the ground motion fields for the selected scenario at the POIs listed in `POIs.txt`.

Preview of an example output file:
<pre>
GMF realizations at Site_LAT:43.2846_LON:12.7778 for Scenario_10: [0.17520797, 0.21844997, 0.093965515, 0.27266037, 0.079073295, 0.09725358, 0.08347481, 0.06693749, 0.005907976, 0.060873847]
GMF realizations at Site_LAT:43.1846_LON:12.8778 for Scenario_10: [0.100996606, 0.35003924, 0.24363522, 0.19941418, 0.15757227, 0.1009447, 0.19146584, 0.06460667, 0.03146108, 0.097111605]
GMF realizations at Site_LAT:43.0846_LON:13.4778 for Scenario_10: [0.18333985, 0.11954803, 0.2914887, 0.050770156, 0.07628956, 0.17871241, 0.10297835, 0.15162756, 0.020328628, 0.04087482]
</pre>

***************************************

**PROB_TOOLS**

`ProbShakemaps` comes with three different tools to generate Probabilistic Shakemaps: 'GetStatistics', 'GetDistributions' and 'EnsemblePlot'. Probabilistic Shakemaps represent different products for visualizing and summarizing the ground-motion predictive distribution at the POIs.

**TOOL: 'GetStatistics'**

* Calculates and save statistics ('Mean','Median','Percentile 10','Percentile 20','Percentile 80','Percentile 90'). If the user does not provide a file with scenarios weights, the scenarios are considered equally probable.
* Plots the calculated statistics at all the selected POIs.

```bash
python ProbShakemaps.py --imt PGA --prob_tool GetStatistics --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 1
```

OUTPUT

* npy files with statistics saved in the `npyFiles` folder
* Statistics map distributions saved in the `STATISTICS` folder

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/STATISTICS/summary_stats_forReadMe.png" alt="SummaryStats" width="90%" height="90%">
</p>

**TOOL: 'GetDistributions'**

Plots the cumulative distribution of the predicted ground-motion values and main statistics at a specific POI together with the ground-motion value recorded at the closest station.

```bash
python ProbShakemaps.py --imt PGA --prob_tool GetDistributions --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 10 --station_file stationlist.json
```

OUTPUT

* `Distr_POI-{POI_idx}.pdf`: Plot of Datum-Ensemble comparison at a given POI
* `POIs_Map.pdf`: Spatial map of the POIs

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/DISTRIBUTIONS/summary_stats_forReadMe.png" alt="DatumEnsemble" width="80%" height="80%">
</p>


**TOOL: 'EnsemblePlot'**

Plots and summarizes the key statistical features of the distribution of predicted ground-motion values at the selected POIs.

```bash
python ProbShakemaps.py --imt PGA --prob_tool EnsemblePlot --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10
```
OUTPUT

* `POIs_Map.pdf`: Spatial map of the POIs
* `Ensemble_Spread_Plot_{imt}.pdf`: Boxplot

<p align="center">
    <img src="https://github.com/angystallone/ProbShakemaps/blob/main/OUTPUT_REPO/summary_stats_forReadMe.png" alt="DatumEnsemble" width="80%" height="80%">
</p>

**POIs SUBSET OPTION**

When using the tools 'QueryHDF5', 'GetStatistics', 'GetDistributions' and 'EnsemblePlot', the user can require to extract a subset of POIs within a maximum distance from the event epicenter following one of these two possible spatial distributions: <ins>random</ins> and <ins>azimuthally uniform</ins>. This changes the command line to):

```bash
python ProbShakemaps.py [...] --pois_subset --n_pois 12 --max_distance 50 --pois_selection_method azimuth_uniform
```
If <ins>azimuthally uniform</ins> is selected, POIs are chosen within a ring in the range ```max_distance +- max_distance/10```.

**MULTIPLE TOOLS AT THE SAME TIME**

```ProbShakemaps``` can handle multiple tools at the same time. Be aware that, in this case, the same settings will apply (ie,```--imt_min```, ```--imt_max```, ```--pois_subset``` etc.). For example:

```bash
python ProbShakemaps.py --imt PGA --prob_tool GetDistributions EnsemblePlot --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 10 --station_file stationlist.json --pois_subset --n_pois 12 --max_distance 50 --pois_selection_method azimuth_uniform
```

**HPC**

```ProbShakemaps```  can be run on a cluster enjoying parallelization. See an example of bash file to run the code on a HPC cluster at [run_code.bash](https://github.com/angystallone/ProbShakemaps/blob/main/run_code.bash). IMPORTANT: the number set at `--ntasks-per-node` must coincide with `num_processes`.


License
-----------------------------

This project is released under the [MIT License](LICENSE).


Contact
-----------------------------

If you need support write to [angela.stallone@ingv.it](mailto:angela.stallone@ingv.it).
