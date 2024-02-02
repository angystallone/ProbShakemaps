# ProbShakemap

`ProbShakemap` is a Python toolbox that propagates source uncertainty from an ensemble of earthquake scenarios to ground motion predictions at a grid of target points. It accounts for model uncertainty by accommodating multiple GMMs and their inherent variability. The package includes `SeisEnsMan`, a tool for generating the ensemble of event-compatible source scenarios. The output consists of a set of products aiding the user to explore and visualize the predictive distribution of ground motion at each target point. Designed for Urgent Computing applications.

Dependencies
------------

 * [Shakemap](https://github.com/DOI-USGS/ghsc-esi-shakemap)
 * [OpenQuake](https://github.com/gem/oq-engine/blob/master/README.md)
 * [Docker](https://www.docker.com/)
 
Command line usage
------------------
<pre>
usage: ProbShakemap.py [-h] [--imt {PGA,PGV,SA(0.3),SA(1.0),SA(3.0}]
                       [--tool {StationRecords,Save_Output,QueryHDF5}]
                       [--prob_tool {GetStatistics,GetDistributions,EnsemblePlot} [{GetStatistics,GetDistributions,EnsemblePlot} ...]]
                       [--numGMPEsRealizations NUMGMPESREALIZATIONS] [--num_processes NUM_PROCESSES]
                       [--imt_min IMT_MIN] [--imt_max IMT_MAX] [--station_file STATION_FILE]
                       [--scenario SCENARIO] [--pois_file POIS_FILE] [--deg_round DEG_ROUND]
                       [--pois_subset] [--n_pois N_POIS] [--max_distance MAX_DISTANCE]
                       [--pois_selection_method {random,azimuth_uniform}]
                       [--fileScenariosWeights FILESCENARIOSWEIGHTS]

ProbShakemap Toolbox

optional arguments:
  -h, --help            show this help message and exit

input params:
  --imt {PGA,PGV,SA(0.3),SA(1.0),SA(3.0)}
                        Intensity measure type (IMT)
  --tool {StationRecords,Save_Output,QueryHDF5}
                        Tool(s) to use
  --prob_tool {GetStatistics,GetDistributions,EnsemblePlot} [{GetStatistics,GetDistributions,EnsemblePlot} ...]
                        ProbShakemap Tool(s) to use
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
---------------

1) <ins>INGV shakemap Docker Image</ins> --> [INGV shakemap Docker Image](https://github.com/INGV/shakemap/releases/tag/v4.1.3) is the INGV configuration of the [USGS Shakemap Docker Image](https://github.com/DOI-USGS/ghsc-esi-shakemap) incorporating specific GMMs and Vs30 map for the Italian region. Except for that, the two products are equivalent. See below for instructions on how to install the INGV shakemap Docker Image.
2) <ins>POIs file</ins> --> two space-separated columns .txt file with LAT and LON of the POIs. The file must be put in the folder `INPUT_FILES`. 
3) <ins>input_file.txt</ins> --> file containing the inputs required by `OpenQuake` and `Shakemap`. The file must be put in the folder `INPUT_FILES` (do not rename it). Be sure to set ID_Event equal to the event_id folder name (see Setting ProbShakemap section below).

> * TectonicRegionType: as defined in OpenQuake tectonic regionalisation.
> * Magnitude_Scaling_Relationship: as required from openquake.hazardlib.scalerel.
> * Rupture_aratio: rupture aspect ratio as required from openquake.hazardlib.geo.surface.PlanarSurface.from_hypocenter
> * ID_Event: Shakemap ID of the event.
> * Vs30file: GMT .grd Vs30 file; if not provided, set it to None. Default Vs30 value (760 m/s) will be used instead.
> * CorrelationModel: as required from openquake.hazardlib.correlation.
> * CrosscorrModel: as required from openquake.hazardlib.cross_orrelation.
> * vs30_clustering: `True` value means that Vs30 values are expected to show clustering (as required from openquake.hazardlib.correlation).
> * truncation_level: number of standard deviations for truncation of the cross-correlation model distribution (as required from openquake.hazardlib.cross_correlation).

4) <ins>fileScenariosWeights.txt</ins> --> File with scenarios weights (optional).  The file must be put in the folder `INPUT_FILES` (do not rename it). 


INSTALLATION
------------

**Set ProbShakemap**

Clone the INGV shakemap GitHub repository (tag v4.1.3) into your working directory:

```bash
git clone --branch v4.1.3 https://github.com/INGV/shakemap.git
```

The folder `shakemap/data/shakemap_profiles/world/data` includes, as an example, the event-id folder for Norcia earthquake (`8863681/current`). The event-id folder contains the file `event.xml`, with basic information about the event. 
You need to create an `event-id/current` folder for each new event and provide the corresponding `event.xml` file. The latter can be built easily: start from the `event.xml` file provided for the Norcia example and then edit latitude, longitude, magnitude and time, the only information needed by `SeisEnsMan` to download the event QUAKEML file (see below). Make sure the event-id is the same you provided in `input_file.txt`. The Vs30 file (`global_italy_vs30_clobber.grd`) is placed in the Docker folder `/home/shake/shakemap_data/vs30` after building the image. The file includes a specific Vs30 model for italy (Michelini et al., 2020). 

Start Docker (download it from [here](https://www.docker.com/)) and build the shakemap Docker Image:

```bash
cd shakemap
DOCKER_BUILDKIT=1 docker build --no-cache --build-arg ENV_UID=$(id -u) --build-arg ENV_GID=$(id -g) --tag shakemap4:4.1.3 .
```

Download `ProbShakemap`:

```bash
git clone https://github.com/INGV/ProbShakemap.git
```
The `ProbShakemap` folder contains the input file, the list of scenarios and the POIs file related to the Amatrice earthquake example.
Move the folder content to the 'world' folder (needed to preserve all the files after shutting down Docker):

```bash
mv ./ProbShakemap/* ./data/shakemap_profiles/world/ && rm -rf ./ProbShakemap
```

**Install SeisEnsMan**

`SeisEnsMan` generates an ensemble of N earthquake source scenarios that are compatible with the event under consideration, given the past seismicity in the region and the known faults. It utilizes the information from the `event.xml` file to automatically download the event's QUAKEML file and generate the `event_stat.json` file. The latter contains all the necessary parameters for creating the ensemble of scenarios. Examples of event-specific JSON files can be found in the `SeisEnsManV2/IO/EarlyEst` folder. 

To install all Python libraries required by `SeisEnsMan`, first create and activate the environment SeisEnsMan:

```bash
python -m venv SeisEnsMan
```

On macOS and Linux:
```bash
source SeisEnsMan/bin/activate
```

On Windows:
```bash
SeisEnsMan\Scripts\activate
```

Then use the file `requirements.txt` provided in the folder `SeisEnsManV2` to install the required libraries:

```bash
python3 -m pip install -r requirements.txt
```

HOW TO RUN
----------

**Generate the scenarios ensemble**

Create the `event-id/current` folder for the event and provide the corresponding `event.xml` file. This will be used by `SeisEnsMan` to download the event QUAKEML file needed for generating the ensemble of event-compatible scenarios.
Activate the environment SeisEnsMan and move to `SeisEnsManV2` directory in `path/to/shakemap/data/shakemap_profiles/world/`. Then run the following command (set the `--nb_scen` parameter to the desired number of scenarios in the ensemble): 

```bash
./line_call.sh
```

After being generated, the ensemble of scenarios is saved in `INPUT_FILES/ENSEMBLE` folder, ready to be queried by `ProbShakemap`. Any other old file has been moved to the `BACKUP` folder.
Before running `ProbShakemap`, make sure to deactivate the environment SeisEnsMan:

```bash
deactivate
```

**Run ProbShakemap**

Start Docker and move back to `shakemap` directory, then run:

```bash
docker run -it --rm -v $(pwd)/data/shakemap_profiles:/home/shake/shakemap_profiles -v $(pwd)/data/local:/home/shake/.local --entrypoint=bash shakemap4:4.1.3
sm_profile -l
cd /home/shake/shakemap_profiles/world
```

`ProbShakemap` comes with three utility tools: `StationRecords`, `Save_Output` and `QueryHDF5`. 

**TOOL: StationRecords**

Plot data from `Shakemap` file `stationlist.json` (the file must be placed in the `event-id/current` folder). 

```bash
python ProbShakemap.py --imt PGA --tool StationRecords --imt_min 0.01 --imt_max 10 --station_file stationlist.json
```
OUTPUT

`Data_stationfile_{imt}.pdf`: Plot data from Shakemap .json station file for the selected IMT (PGA in the example).

<p align="center">
    <img src="https://github.com/INGV/ProbShakemap/blob/main/OUTPUT_REPO/Data_stationfile_PGA_2.png" alt="Data_stationfile_PGA" width="60%" height="60%">
</p>


**TOOL: Save_Output**

Run the probabilistic analysis and save the output to a .HDF5 file with the following hierarchical structure.

scenario --> POI --> GMPEs realizations

```bash
python ProbShakemap.py --imt PGA --tool Save_Output --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10
```

OUTPUT

`SIZE_{num_scenarios}_ENSEMBLE_{IMT}.hdf5`


**TOOL: QueryHDF5**

Navigate and query the .HDF5 file.

```bash
python ProbShakemap.py --tool QueryHDF5 --imt PGA --scenario 10 --pois_file POIs.txt
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

`ProbShakemap` comes with three 'prob tools': `GetStatistics`, `GetDistributions` and `EnsemblePlot`. The outputs are intended to assist the user in exploring the ground-motion predictive distributions at a set of POIs.

**TOOL: GetStatistics**

Calculate and save statistics of the predictive distribution. Plot the calculated statistics at the POIs.

```bash
python ProbShakemap.py --imt PGA --prob_tool GetStatistics --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 1
```

OUTPUT

* npy files with statistics saved in the `npyFiles` folder
* map distributions of statistics in `vector_stat.npy` saved in the `STATISTICS` folder

The `npyFiles` folder contains:
* `vector.npy`: a 2D array that stores the ground motion distributions across all POIs. The array has dimensions (`num_pois`, `num_GMPEsRealizations` * `num_scenarios`), where `num_GMPEsRealizations` represents the number of realizations per scenario, and `num_scenarios` is the total number of scenarios in the ensemble; 
* `thresholds_distrib.npy`: a 2D array representing the probabilistic distributions of ground motion across POIs. The array has dimensions (`num_pois`, 6,000), where 6,000 is the number of intervals into which the range of ground motion values has been discretized. Probabilities within each interval are weighted by the scenario weights and then aggregated across all scenarios in the ensemble;
* `thresholds_stat.npy`: dictionary of statistics derived from the distributions in `thresholds_distrib.npy`. For each Point of Interest (POI), it returns: 'Weighted Mean', 'Arithmetic Mean', 'Median','Percentile 10','Percentile 20','Percentile 80','Percentile 90'; 
* `vector_stat.npy`: dictionary of weighted statistics computed from the distributions in `vector.npy`. For each Point of Interest (POI), it returns: 'Mean', 'Median','Percentile 10','Percentile 20','Percentile 80','Percentile 90','Percentile 5','Percentile 95','Percentile 2.5','Percentile 97.5';
* `weight.npy`: a 2D array that stores the normalized weights across all POIs. The array dimensions match those of `vector.npy`, which are (`num_pois`, `num_GMPEsRealizations` * `num_scenarios`).


<p align="center">
    <img src="https://github.com/INGV/ProbShakemap/blob/main/OUTPUT_REPO/STATISTICS/summary_stats_forReadMe.png" alt="SummaryStats" width="90%" height="90%">
</p>

**TOOL: GetDistributions**

Plot the cumulative distribution of the predicted ground-motion values and main statistics at a specific POI together with the ground-motion value recorded at the closest station (or at a POI coincident with the station, if available).

Note: the `Shakemap` file `stationlist.json` must be placed in the `event-id/current` folder. 

```bash
python ProbShakemap.py --imt PGA --prob_tool GetDistributions --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 10 --station_file stationlist.json
```

OUTPUT

* `POIs_Map.pdf`: Spatial map of the POIs
* `Distr_POI-{POI_idx}.pdf`: Plot of Datum-Ensemble comparison at a given POI

<p align="center">
    <img src="https://github.com/INGV/ProbShakemap/blob/main/OUTPUT_REPO/POIs_Map.png" alt="DatumEnsemble" width="25%" height="25%">
</p>

<p align="center">
    <img src="https://github.com/INGV/ProbShakemap/blob/main/OUTPUT_REPO/DISTRIBUTIONS/summary_stats_forReadMe.png" alt="DatumEnsemble" width="90%" height="90%">
</p>


**TOOL: EnsemblePlot**

Plot and summarize the key statistical features of the distribution of predicted ground-motion values at the POIs.

```bash
python ProbShakemap.py --imt PGA --prob_tool EnsemblePlot --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10
```

OUTPUT

* `POIs_Map.pdf`: Spatial map of the POIs
* `Ensemble_Spread_Plot_{imt}.pdf`: Boxplot

<p align="center">
    <img src="https://github.com/INGV/ProbShakemap/blob/main/OUTPUT_REPO/Ensemble_Spread_Plot.png" alt="DatumEnsemble" width="50%" height="50%">
</p>

**POIs SUBSET OPTION**

When using the tools `QueryHDF5`, `GetStatistics`, `GetDistributions` and `EnsemblePlot`, you can require to extract a subset of POIs within a maximum distance from the event epicenter following one of the following spatial distributions: <ins>random</ins> and <ins>azimuthally uniform</ins>. This changes the command line to:

```bash
python ProbShakemap.py [...] --pois_subset --n_pois 12 --max_distance 50 --pois_selection_method azimuth_uniform
```
If <ins>azimuthally uniform</ins> is selected, POIs are chosen within a ring in the range `max_distance +- max_distance/10`.

**MULTIPLE TOOLS AT THE SAME TIME**

`ProbShakemap` can handle multiple tools at the same time. Be aware that, in this case, the same settings will apply (ie,`--imt_min`, `--imt_max`, `--pois_subset` etc.).

```bash
python ProbShakemap.py --imt PGA --prob_tool GetDistributions EnsemblePlot --num_processes 8 --pois_file POIs.txt --numGMPEsRealizations 10 --imt_min 0.001 --imt_max 10 --station_file stationlist.json --pois_subset --n_pois 12 --max_distance 50 --pois_selection_method azimuth_uniform
```

**HPC**

`ProbShakemap`  can be run on a cluster enjoying parallelization. See an example of bash file to run the code on a HPC cluster at [run_code.bash](https://github.com/angystallone/ProbShakemap/blob/main/run_code.bash). IMPORTANT: the number set at `--ntasks-per-node` must coincide with `num_processes`.


License
-------

This project is released under the [MIT License](LICENSE).


Contact
--------

If you need support write to [angela.stallone@ingv.it](mailto:angela.stallone@ingv.it).


Contributions & Acknowledgements
--------------------------------

Jacopo Selva coded the `GetStatistics` tool; Louise Cordrie authored the `SeisEnsMan` tool and tested `ProbShakemap` on the INGV-Bologna ADA cluster.
I thank Valentino Lauciani for testing and developing the INGV Shakemap Docker and Licia Faenza for testing ProbShakemap. 

To cite ProbShakemap
--------------------
