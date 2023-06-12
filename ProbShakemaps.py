import tools
import os
import numpy
import argparse
import time

start_time = time.time()

def run_main():
    
    # set the random seed for reproducibility
    seed = 0
    numpy.random.seed(seed)

    ProbCalc = tools.Main(args.imt, args.pois_file,
                        args.numGMPEsRealizations, args.num_processes)
    
    prob_output = ProbCalc.run_prob_analysis()

    return prob_output


# Define command-line arguments
parser = argparse.ArgumentParser(description='ProbShakemap Toolbox')
input_params = parser.add_argument_group('input params')
input_params.add_argument('--imt', choices=['PGA', 'PGV'], help='Intensity measure type (IMT)')
input_params.add_argument('--tool', choices=['StationRecords', 'Save_Output', 'QueryHDF5'], help='Tool(s) to use')
input_params.add_argument('--prob_tool', choices=['GetStatistics', 'GetDistributions', 'EnsemblePlot'], nargs='+', help='ProbShakemaps Tool(s) to use')
input_params.add_argument('--numGMPEsRealizations', type=int, help='Total number of GMPEs random samples')
input_params.add_argument('--num_processes', type=int, default=1, help='Number of CPU cores for code parallelization')
input_params.add_argument('--imt_min', type=float, help='Minimum value for the selected IMT (for plot only)')
input_params.add_argument('--imt_max', type=float, help='Maximum value for the selected IMT (for plot only)')
input_params.add_argument('--station_file', help='Shakemap .json station file')
input_params.add_argument('--scenario', type=float, help='Scenario number')
input_params.add_argument('--pois_file', help='Filename with latitude and longitude of POIs')
input_params.add_argument('--deg-round', type=float, default=5, help='Rounding precision for latitude and longitude')
input_params.add_argument('--pois_subset', action='store_true', default=False, help='Extract a subset of POIs')
input_params.add_argument('--n_pois', type=int, default=10, help='Number of POIs in the subset')
input_params.add_argument('--max_distance', type=int, default=200, help='Max distance from epicenter of POIs in the subset')
input_params.add_argument('--pois_selection_method', choices=['random', 'azimuth_uniform'], help='Selection method for the POIs of the subset')
input_params.add_argument('--fileScenariosWeights', default="", help='File with scenarios weights')

# Parse command-line arguments
args = parser.parse_args()

# Get params
params = tools.get_params()
Lat_Event = float(params["Lat_Event"])
Lon_Event = float(params["Lon_Event"])
event_dir = params["event_dir"]
EnsembleSize = params['Ensemble_Size']

#################################
##### CHECK STATION RECORDS #####
#################################

if args.tool == 'StationRecords':

    if args.imt_max is None:
        raise TypeError("Missing required argument 'imt_max'")
    if args.imt_min is None:
        raise TypeError("Missing required argument 'imt_min'")
    if args.station_file is None:
        raise TypeError("Station .json file from shakemap not available")

    StationRecords = tools.StationRecords(args.imt, args.imt_min, args.imt_max, 
                                            args.deg_round, args.station_file)
    StationRecords.plot()


#################################
##### WRITE OUTPUT FILE #########
#################################

elif args.tool == 'Save_Output':

    if args.imt is None:
        raise TypeError("Missing required argument 'imt'")
    if args.pois_file is None:
                    raise TypeError("Missing required argument 'pois_file'")
    if args.numGMPEsRealizations is None:
        raise TypeError("Missing required argument 'numGMPEsRealizations'")
        
    
    prob_output = run_main()

    SiteGmf = prob_output["SiteGmf"]
    keys_scen = prob_output["keys_scen"]
    keys_sites = prob_output["keys_sites"]

    Save_Output = tools.Write(args.imt, EnsembleSize, keys_scen, 
                              SiteGmf, keys_sites, args.num_processes)
    Save_Output.write_output()


#################################
##### QUERY OUTPUT FILE #########
#################################

elif args.tool == 'QueryHDF5':

    outfile_dir = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/")
    outputfile = [name for name in os.listdir(outfile_dir) if name != ".DS_Store"][0]
    if not outputfile:
        print("WARNING: No output file found --> run 'Save_Output' first!")
    else:
        if args.scenario is None:
            raise TypeError("Missing required argument 'scenario'")

    QueryHDF5 = tools.QueryHDF5(args.scenario, args.pois_file,
                                args.pois_subset, args.n_pois, args.max_distance, 
                                args.pois_selection_method, Lon_Event, Lat_Event)
    QueryHDF5.print_info()


#################################
###### PROB SHAKEMAPS TOOLS #####
#################################

if args.prob_tool:

    run_main_flag = False  # Flag variable to track if Main() has been executed
    
    for tool in args.prob_tool:
            
        if tool == 'GetStatistics':
            
            if args.imt_max is None:
                raise TypeError("Missing required argument 'imt_max'")
            if args.imt_min is None:
                raise TypeError("Missing required argument 'imt_min'")
            if args.pois_subset and not args.pois_selection_method:
                raise TypeError("Missing required argument 'pois_selection_method'")
            if args.pois_selection_method == 'azimuth_uniform':
                if args.n_pois % 4 > 0:
                    raise TypeError("Select a number of POIs divisible by 4")

            if run_main_flag == False:

                if args.imt is None:
                    raise TypeError("Missing required argument 'imt'")
                if args.pois_file is None:
                    raise TypeError("Missing required argument 'pois_file'")
                if args.numGMPEsRealizations is None:
                    raise TypeError("Missing required argument 'numGMPEsRealizations'")
                    
                prob_output = run_main()
                SiteGmf = prob_output["SiteGmf"]

                run_main_flag = True    
                
            GetStatistics = tools.GetStatistics(SiteGmf, EnsembleSize, Lon_Event, Lat_Event, args.numGMPEsRealizations, event_dir, args.imt, 
                                                args.imt_min, args.imt_max, args.fileScenariosWeights, 
                                                args.pois_file, args.pois_subset, args.n_pois, args.max_distance, 
                                                args.pois_selection_method, args.deg_round)

            GetStatistics.save_statistics()
            GetStatistics.plot_statistics()


        elif tool == 'GetDistributions':

            if args.imt_max is None:
                raise TypeError("Missing required argument 'imt_max'")
            if args.imt_min is None:
                raise TypeError("Missing required argument 'imt_min'")
            if args.station_file is None:
                raise TypeError("Station .json file from shakemap not available")
            if args.pois_subset and not args.pois_selection_method:
                raise TypeError("Missing required argument 'pois_selection_method'")
            if args.pois_selection_method == 'azimuth_uniform':
                if args.n_pois % 4 > 0:
                    raise TypeError("Select a number of POIs divisible by 4")

            if run_main_flag == False:

                if args.imt is None:
                    raise TypeError("Missing required argument 'imt'")
                if args.pois_file is None:
                    raise TypeError("Missing required argument 'pois_file'")
                if args.numGMPEsRealizations is None:
                    raise TypeError("Missing required argument 'numGMPEsRealizations'")
                    
                prob_output = run_main()
                SiteGmf = prob_output["SiteGmf"]

                run_main_flag = True    

            GetDistributions = tools.GetDistributions(SiteGmf, EnsembleSize, Lon_Event, Lat_Event, args.numGMPEsRealizations, event_dir, args.imt, args.station_file,
                                                    args.imt_min, args.imt_max, args.fileScenariosWeights, 
                                                    args.pois_file, args.pois_subset, args.n_pois, args.max_distance, 
                                                    args.pois_selection_method, args.deg_round)

            GetDistributions.plot_distributions() 


        elif tool == 'EnsemblePlot':

            if args.pois_subset and not args.pois_selection_method:
                raise TypeError("Missing required argument 'pois_selection_method'")    
            if args.pois_selection_method == 'azimuth_uniform':
                if args.n_pois % 4 > 0:
                    raise TypeError("Select a number of POIs divisible by 4")

            if run_main_flag == False:

                if args.pois_file is None:
                    raise TypeError("Missing required argument 'pois_file'")
                if args.imt is None:
                    raise TypeError("Missing required argument 'imt'")
                if args.numGMPEsRealizations is None:
                    raise TypeError("Missing required argument 'numGMPEsRealizations'")
                    
                prob_output = run_main()
                SiteGmf = prob_output["SiteGmf"]
                run_main_flag = True    

            EnsemblePlot = tools.EnsemblePlot(SiteGmf, args.imt, Lon_Event, Lat_Event, EnsembleSize, args.pois_file, args.pois_subset, 
                                                args.n_pois, args.max_distance, args.deg_round, args.pois_selection_method)

            EnsemblePlot.plot()

        else:

            print('No tool specified')

print("********* DONE! *******")
print("--- %s seconds ---" % (time.time() - start_time))




