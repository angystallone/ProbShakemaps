import tools
import main
import json
import os
import numpy
import argparse

params_file = 'params.json'

# Define command-line arguments
parser = argparse.ArgumentParser(description='ProbShakemap Toolbox')
input_params = parser.add_argument_group('input params')
input_params.add_argument('--task', choices=['RunProbAnalysis', 'GenerateProbShakemap'], required=True, help='Task to perform')
input_params.add_argument('--imt', choices=['PGA', 'PGV'], help='Intensity measure type (IMT)')
input_params.add_argument('--tool', choices=['StationRecords', 'QueryHDF5', 'GetStatistics', 'EnsemblePlot'], help='Tool to use')
input_params.add_argument('--imt_min', type=float, help='Minimum value for the selected IMT (for plot only)')
input_params.add_argument('--imt_max', type=float, help='Maximum value for the selected IMT (for plot only)')
input_params.add_argument('--station_file', help='Shakemap .json station file')
input_params.add_argument('--scenario', type=float, help='Scenario number')
input_params.add_argument('--pois_file', help='Filename with latitude and longitude of POIs')
input_params.add_argument('--ensemble_number', type=int, default=1, help='Ensemble sample number')
input_params.add_argument('--deg-round', type=float, default=5, help='Rounding precision for latitude and longitude')
input_params.add_argument('--pois_subset', action='store_true', default=False, help='Extract a subset of POIs')
input_params.add_argument('--n_pois', type=int, default=10, help='Number of POIs in the subset')
input_params.add_argument('--max_distance', type=int, default=200, help='Max distance from epicenter of POIs in the subset')
input_params.add_argument('--pois_selection_method', choices=['random', 'azimuth_uniform'], help='Selection method for the POIs of the subset')
input_params.add_argument('--fileScenariosWeights', default="", help='File with scenarios weights')

# Parse command-line arguments
args = parser.parse_args()

if args.task is None:
    raise TypeError("Missing required argument 'task'")

# Run probabilistic analysis if requested
if args.task == 'RunProbAnalysis':
    
    print("***** Running probabilistic analysis *****")

    filename = os.path.join(os.getcwd(), "OUTPUT/", params_file)
    # Run main.run_prob_analysis() and save the params to the file
   
    # set the random seed for reproducibility
    seed = 0
    numpy.random.seed(seed)

    params = main.run_prob_analysis()
    with open(filename, 'w') as f:
        json.dump(params, f)

# Generate prob shakemaps
if args.task == 'GenerateProbShakemap':
    if args.imt is None:
        raise TypeError("Missing required argument 'imt'")
          
    print("IMT = ", args.imt)
    print("Ensemble sample number = ", args.ensemble_number)     
    print("***** Loading params file *****")
    filename = os.path.join(os.getcwd(), "OUTPUT/", params_file)
    with open(filename, 'r') as f:
        params = json.load(f)

    Lat_Event = float(params[args.imt][str(args.ensemble_number - 1)]['Lat_Event'])
    Lon_Event = float(params[args.imt][str(args.ensemble_number - 1)]['Lon_Event'])
    event_dir = params[args.imt][str(args.ensemble_number - 1)]['event_dir'] 
    NumGMPEsRealizations = int(params[args.imt][str(args.ensemble_number - 1)]['NumGMPEsRealizations'])
    output = params[args.imt][str(args.ensemble_number - 1)]['HDF5_Filename']
    print("Lat_Event = ", Lat_Event)
    print("Lon_Event = ", Lon_Event)
    print("NumGMPEsRealizations = ", NumGMPEsRealizations)


    # Run selected tool
    if args.tool == 'StationRecords':

        if args.imt is None:
            raise TypeError("Missing required argument 'imt'")
        if args.imt_max is None:
            raise TypeError("Missing required argument 'imt_max'")
        if args.imt_min is None:
            raise TypeError("Missing required argument 'imt_min'")
        if args.station_file is None:
            raise TypeError("Station .json file from shakemap not available")

        StationRecords = tools.StationRecords(event_dir, args.imt, args.imt_min, args.imt_max, 
                                              args.deg_round, args.station_file)
        StationRecords.plot()
        
    elif args.tool == 'QueryHDF5':

        if args.scenario is None:
            raise TypeError("Missing required argument 'scenario'")
        if args.pois_file is None:
            raise TypeError("Missing required argument 'pois_file'")

        QueryHDF5 = tools.QueryHDF5(output, Lon_Event, Lat_Event, args.scenario, args.pois_file, 
                                    args.pois_subset, args.n_pois, args.max_distance, args.pois_selection_method)
        QueryHDF5.print_info()

    elif args.tool == 'GetStatistics':

        if args.imt is None:
            raise TypeError("Missing required argument 'imt'")
        if args.imt_max is None:
            raise TypeError("Missing required argument 'imt_max'")
        if args.imt_min is None:
            raise TypeError("Missing required argument 'imt_min'")
        if args.station_file is None:
            raise TypeError("Station .json file from shakemap not available")
        if args.pois_file is None:
            raise TypeError("Missing required argument 'pois_file'")
   
        GetStatistics = tools.GetStatistics(output, Lon_Event, Lat_Event, event_dir, args.imt, args.station_file,
                                             args.imt_min, args.imt_max, args.fileScenariosWeights, 
                                             args.pois_file, args.pois_subset, args.n_pois, args.max_distance, 
                                             args.pois_selection_method, args.deg_round)
    
        #GetStatistics.save_statistics()
        GetStatistics.plot_statistics()

    elif args.tool == 'GetDistributions':

        if args.imt is None:
            raise TypeError("Missing required argument 'imt'")
        if args.imt_max is None:
            raise TypeError("Missing required argument 'imt_max'")
        if args.imt_min is None:
            raise TypeError("Missing required argument 'imt_min'")
        if args.station_file is None:
            raise TypeError("Station .json file from shakemap not available")
        if args.pois_file is None:
            raise TypeError("Missing required argument 'pois_file'")
   
        GetDistributions = tools.GetDistributions(output, Lon_Event, Lat_Event, event_dir, args.imt, args.station_file,
                                             args.imt_min, args.imt_max, args.fileScenariosWeights, 
                                             args.pois_file, args.pois_subset, args.n_pois, args.max_distance, 
                                             args.pois_selection_method, args.deg_round)
    
        GetDistributions.plot_distributions() 

    elif args.tool == 'EnsemblePlot':

        if args.imt is None:
            raise TypeError("Missing required argument 'imt'")
        if args.pois_file is None:
            raise TypeError("Missing required argument 'pois_file'")

        EnsemblePlot = tools.EnsemblePlot(output, args.imt, Lon_Event, Lat_Event, args.pois_file, args.pois_subset, 
                                          args.n_pois, args.max_distance, args.deg_round, args.pois_selection_method)

        EnsemblePlot.plot()

    else:

        print('No tool specified')
