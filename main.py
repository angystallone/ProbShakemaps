# %%

# -------------------------------- RUN THE SHAKEMAP CONTAINER AND THE CODE WITHIN IT -------------------------------- #

# cd shakemap4/
# docker run -it --rm -t -p 8888:8888 -v $(pwd)/data/shakemap_profiles:/home/shake/shakemap_profiles -v $(pwd)/data/shakemap_data:/home/shake/shakemap_data -v $(pwd)/data/local:/home/shake/.local --entrypoint=bash shakemap4
# sm_profile -l (--> Profile: world **Current Profile**, Install Path: /home/shake/shakemap_profiles/world/install, Data Path: /home/shake/shakemap_profiles/world/data)
# cd /home/shake/shakemap_profiles/world
# python GMPEs.py

# NOTE: Since the Vs30 file was already downloaded, I had to change its relative path in model.conf

# %%

from distutils.command import install
from openquake.hazardlib.source import BaseRupture
from openquake.hazardlib.geo import Point, surface
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.gsim.base import ContextMaker
from openquake.hazardlib.calc.gmf import GmfComputer
from openquake.hazardlib import valid

import shakemap
from shakemap.utils.config import get_config_paths
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.select import SelectModule
from shakemap.coremods.model import ModelModule
from mapio.gmt import GMTGrid

from configobj import ConfigObj

import numpy
import os
import importlib
import time

import xml.etree.ElementTree as ET
import h5py

from matplotlib import pyplot as plt
from matplotlib.colorbar import cm

import multiprocessing
from multiprocessing import Pool, Process, Lock, Manager
import itertools

from functools import partial
import config

#################################################################


# To sort filenames by sample number
def grab_sample_number(x):
    return(x.split('.')[-2:])


def process_scenario(scen, Ensemble_Scenarios, imts, msr, rupture_aratio, 
                     tectonicRegionType, context_maker, Site_Collection, 
                     correlation_model, crosscorr_model, gmpes, Weighted_Num_Realiz):

    # Get scenario index 
    k = Ensemble_Scenarios.index(scen)

    # Initialize GMF output for this scenario
    gmf = [[] for _ in range(len(imts))]

    Mag = float(scen[0])
    Hypocenter = Point(float(scen[1]), float((scen[2])), float(scen[3]))
    Rake = float(scen[6])
    Strike = float(scen[4])
    Dip = float(scen[5])

    planar_surface = surface.PlanarSurface.from_hypocenter(
    hypoc=Hypocenter,
    msr=msr(), 
    mag=Mag,
    aratio=rupture_aratio,
    strike=Strike,
    dip=Dip,
    rake=Rake,
    )

    source = BaseRupture(
    mag=Mag,
    rake=Rake,
    tectonic_region_type=tectonicRegionType,
    hypocenter=Hypocenter,
    surface=planar_surface
    )

    ctx = context_maker.get_ctxs([source], Site_Collection)
    gc = GmfComputer(
        source, 
        Site_Collection, 
        context_maker,
        correlation_model=correlation_model,
        cross_correl=crosscorr_model
        ) 

    #print(gc.seed)

    mean_and_stdev = context_maker.get_mean_stds(ctx)  # Get an array of shape (4, G, M, N) with mean and stddevs 
        
    # Loop over GMPEs    
    for g, gmpe in enumerate(gmpes):
        # 'gc.compute' --> Compute gmf and returns an array of shape (num_imts, num_sites, num_events) with the sampled gmf, and two arrays with shape 
        # (num_imts, num_events): sig for tau and eps for the random part
        gf = gc.compute(gmpe, Weighted_Num_Realiz[g], mean_and_stdev[:, g, :, :])  

        # Take only first array output from gc.compute (shape: (num_imts, num_sites, num_events)) and loop over IMTS 
        for m in range(len(imts)):

            gmf[m].append(gf[0][m].view())

    return k, gmf


def write_hdf5_file(m, j, i, chunk_size, num_processes, EnsembleSize, subdir1, imts, keys_scen, sites, gmpes, 
                    GMPEsRealizationsForProbShakeMap_AllGMPEs, SiteGmf, keys_sites, lock):

    print(f"Writing chunk {i + 1} ...")
    
    k_start = i * chunk_size
    k_end = (i+1) * chunk_size
    # adjust k_start and k_end for the last chunk
    if i == num_processes - 1:
        k_start = EnsembleSize - chunk_size
        k_end = EnsembleSize 

    filename = subdir1 + "/" + "SIZE_%d" % EnsembleSize + "_ENSEMBLE_%d_%s_chunk%d-%d.hdf5" % (j + 1, imts[m], i, i+1)

    # One .hdf5 file per ensemble and per imts
    with h5py.File(filename,'w') as h5f:
        # Loop over source scenarios within each sampled ensemble
        for k in range(k_start, k_end):
            # One group per scenario
            grp = h5f.create_group(keys_scen[k])

            # Loop over sites
            for s in range(len(sites)):
                SiteGmfGMPE = []
                for g in range(len(gmpes)): 
                    SiteGmfGMPE.append(GMPEsRealizationsForProbShakeMap_AllGMPEs[j][k][m][g][s])

                # Use a lock to synchronize access to the SiteGmf array
                lock.acquire()
                try:
                    SiteGmf[m][j][k][s] = [x for sublist in SiteGmfGMPE for x in sublist]
                finally:
                    lock.release()

                # One dataset per site    
                data_name = keys_sites[s]
                grp.create_dataset(data_name, data=SiteGmf[m][j][k][s], compression="gzip", compression_opts=5)

    print(f"Chunk {i + 1} DONE")


def run_prob_analysis():

    # Load configuration parameters
    config_dict = config.load_config('input_file.txt')

    tectonicRegionType = config_dict['tectonicRegionType']
    mag_scaling = config_dict['mag_scaling']
    rupture_aratio = config_dict['rupture_aratio']
    ID_Event = config_dict['ID_Event']
    POIsfile = config_dict['POIsfile']
    vs30file = config_dict['vs30file']
    CorrelationModel = config_dict['CorrelationModel']
    CrosscorrModel = config_dict['CrosscorrModel']
    vs30_clustering = config_dict['vs30_clustering']
    truncation_level = config_dict['truncation_level']
    NumGMPEsRealizations = config_dict['NumGMPEsRealizations']
    imts = config_dict['imts']
    num_processes = config_dict['num_processes']

    # DEFINE OUT DIR 

    path = os.path.join(os.getcwd(), "OUTPUT")
    if not os.path.exists(path):
        os.makedirs(path)

    subdir1 = path + "/HDF5_FILES" 
    if not os.path.exists(subdir1):
        os.makedirs(subdir1)

    # PRINTING INPUT INFO

    print("TectonicRegionType: " + tectonicRegionType)

    print("Importing " + mag_scaling + " as magnitude scaling relationship")
    module = importlib.import_module('openquake.hazardlib.scalerel')
    msr = getattr(module, mag_scaling)

    print("Event ID: " + ID_Event)

    print("POIs file: " + POIsfile)

    if vs30file == "":
        print("Vs30 file not provided")
    else:
        print("Vs30 file: " + vs30file)    

    print("Importing " + CorrelationModel + " as correlation model")
    module = importlib.import_module('openquake.hazardlib.correlation')
    correlation_model = getattr(module, CorrelationModel)

    print("Importing " + CrosscorrModel + " as crosscorrelation model")
    module = importlib.import_module('openquake.hazardlib.cross_correlation')
    crosscorr_model = getattr(module, CrosscorrModel)

    print("Vs30 clustering: " + str(vs30_clustering))

    print("Truncation level: " + str(truncation_level))

    print("Number of GMPEs realizations per POI: " + str(NumGMPEsRealizations))

    print("Intensity measure types: " + str(imts))

    print("Number of CPU processes: " + str(num_processes))

    # LOAD INFO FROM 'shake_result.hdf'

    # Install dir and event dir
    install_path, data_path = get_config_paths()
    print("Install Path = ", install_path)
    print("Data Path = ", data_path)
    event_dir = os.path.join(data_path, ID_Event, "current")

    if not os.path.exists(event_dir):
        raise NotADirectoryError(f"{event_dir} is not a valid directory.")

    eventxml = os.path.join(event_dir, "event.xml")    
    if not os.path.isfile(eventxml):
        raise FileNotFoundError(f"{eventxml} does not exist.")   

    # Parse the eventxml file
    tree = ET.parse(eventxml)
    root = tree.getroot()
    # Get the latitude and longitude of the event
    Lat_Event = root.attrib['lat']
    Lon_Event = root.attrib['lon']

    # Collects event and configuration data and creates the file shake_data.hdf
    assemble = AssembleModule(ID_Event, comment='Test comment.')
    assemble.execute()

    # Reads the data in shake_data.hdf and produces an interpolated ShakeMap
    model = ModelModule(ID_Event)
    model.execute()

    # Generate model_select.conf file, containing GMPE sets 
    select = SelectModule(ID_Event)
    select.execute()

    # Read the event.xml file and generate a GMPE set for the event based on the event’s residence within, 
    # and proximity to, a set of predefined tectonic regions and user-defined geographic areas

    print("********* RETRIEVING GMPEs *******")

    # Extract GMPE set
    conf_filename = event_dir + '/model_select.conf'
    config_model_select = ConfigObj(conf_filename)
    print("Config filename = ", conf_filename)

    GMPE_Set = []
    for key, _ in config_model_select['gmpe_sets'].items():
        GMPE_Set.append(config_model_select['gmpe_sets'][key]['gmpes'][0])

    print("GMPE Sets selected for this earthquake = ", GMPE_Set)

    # config files
    conf_filename = install_path + '/config/gmpe_sets.conf'
    config_gmpe_sets = ConfigObj(conf_filename)

    conf_filename = install_path + '/config//modules.conf'
    config_modules = ConfigObj(conf_filename)

    # Get GMPEs acronyms
    GMPEs_Acronym, GMPEs_Weights, GMPEs_Names = [], [], []
    for key, _ in config_gmpe_sets['gmpe_sets'].items():
        if key in GMPE_Set:
            GMPEs_Acronym.append(config_gmpe_sets['gmpe_sets'][key]['gmpes'])
            GMPEs_Weights.append(float(config_gmpe_sets['gmpe_sets'][key]['weights']))

    # Get GMPEs from the acronyms
    for key, item in config_modules['gmpe_modules'].items():  
        if key in GMPEs_Acronym:
            print(f"Importing {item[0]}")
            GMPEs_Names.append(item[0])

    # Convert gmpes names into GSIM instances
    gmpes = []
    for elem in GMPEs_Names:
        gmpes.append(valid.gsim(elem))  # OQ equivalent of getattr

    # Load POIs lat and long
    print("********* LOADING POIs *******")

    filename = POIsfile
    path = os.path.join(os.getcwd(), "INPUT_FILES")
    with open(path + '/' + filename, 'r') as f:
        POIs = [x.strip().split() for x in f]

    print("Number of POIs = ", len(POIs))

    # Load LON and LAT of POIs
    LON = numpy.empty(len(POIs))
    LAT = numpy.empty(len(POIs))
    for i in range(len(POIs)):
        LON[i] = float(POIs[i][1])
        LAT[i] = float(POIs[i][0])

    if vs30file:

        print("********* LOADING Vs30 *******")
        # Vs30
        vs30fullname = os.path.join(data_path, 'shakemap_data', 'vs30', vs30file)
        vs30grid = GMTGrid.load(vs30fullname)
        # Interpolate Vs30 values at POIs 
        vs30_POIs = vs30grid.getValue(LAT, LON, method="nearest")

    print("********* DEFINING SITE COLLECTION *******")

    # Define a SiteCollection for all the POIs
    sites = []
    for i in range(len(POIs)):
        site_location = Point(LON[i], LAT[i])
        # If Vs30 file is provided
        if vs30file:
            site = Site(location=site_location, vs30=vs30_POIs[i], vs30measured=False, z1pt0=40., z2pt5=1.0)
        # If Vs30 file is not provided, use default value for Vs30
        else:
            site = Site(location=site_location, vs30=760., z1pt0=40., z2pt5=1.0)
        sites.append(site)
        
    Site_Collection = SiteCollection(sites)
    print(Site_Collection.complete)    

    listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"

    # Get ensemble file list (will contain more than 1 file if multiple ensembles of the same size are available)    
    filelist = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]

    NumSampledEnsembles = len(filelist)
    print("Number of random ensemble samples : " + str(NumSampledEnsembles))

    # Sort filenames by sample number    
    filelist_sorted = []
    for file in sorted(filelist, key = grab_sample_number): 
        filelist_sorted.append(file)

    # Get the number of scenarios
    with open(os.path.join(listscenarios_dir, filelist_sorted[0]), 'r') as f:
        EnsembleSize = 0
        for line in f:
            EnsembleSize += 1

    print("Number of source scenarios to process = ", EnsembleSize)

    # Build contexts

    print("********* BUILDING OPENQUAKE CONTEXTS *******")

    # # Define input parameters for ContextMaker

    imtls = {}
    for imt in imts:
        imtls[imt] = []

    param = dict(imtls=imtls)   

    # Instantiate a ContextMaker object (Note: independent from source and sites!)
    context_maker = ContextMaker(tectonicRegionType, gmpes, param)

    print("********* SAMPLING UNCERTAINTY *******")

    correlation_model = correlation_model(vs30_clustering=vs30_clustering)
    crosscorr_model = crosscorr_model(truncation_level=truncation_level)

    # This is needed to sample proportionally to the weight of the GMPEs
    Weighted_Num_Realiz = []
    for i in range(len(GMPEs_Names)):
        Weighted_Num_Realiz.append(round(NumGMPEsRealizations * GMPEs_Weights[i]))

    # Check
    a = sum(Weighted_Num_Realiz) == NumGMPEsRealizations
    if a == False:
        raise Exception(f"Something wrong with the number of GMPEs realizations.") 

    start_time = time.time()

    # Sample from the total variability of ground motion taking into account both inter- and intra-event variability (for one source scenario only)
    # gmf = exp(mu + crosscorel(tau) + spatialcorrel(phi)) --> See: https://docs.openquake.org/oq-engine/advanced/latest/event_based.html#correlation-of-ground-motion-fields

    GMPEsRealizationsForProbShakeMap_AllGMPEs = numpy.zeros((NumSampledEnsembles, EnsembleSize, len(imts)), dtype=object)

    for j in range(NumSampledEnsembles):

        f = open(os.path.join(listscenarios_dir, filelist_sorted[j]), 'r')
        print("Processing file:", filelist_sorted[j])
        Ensemble_Scenarios = []
        for k, line in enumerate(f):
            scen = line.strip().split(' ')
            Ensemble_Scenarios.append(scen)  

        # Set the random seed for reproducibility in OpenQuake GmfComputer
        seed = 0
        numpy.random.seed(seed)

        ################################
        # SETTING MULTIPROCESSING PARAMS
        ################################

        # Number of CPU cores to use
        num_processes = num_processes
        chunk_size_default = int(EnsembleSize/num_processes) # size of each chunk of scenarios
        print("Chunk size = ", chunk_size_default)
        last_chunk_size = chunk_size_default + EnsembleSize - num_processes * chunk_size_default # size of the last chunk
        print("Last chunk size = ", last_chunk_size)

        # Create pool of worker processes
        with Pool(processes=num_processes) as pool:
            results = []

            # iterate over processes
            for i in range(num_processes):
                if i == num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default

                start_idx = i * chunk_size
                end_idx = (i+1) * chunk_size
                # adjust k_start and k_end for the last chunk
                if i == num_processes - 1:
                    start_idx = EnsembleSize - chunk_size
                    end_idx = EnsembleSize 

                chunk = Ensemble_Scenarios[start_idx:end_idx]

                chunk_results = []
                for scenario in chunk:
                    result = process_scenario(scen=scenario, Ensemble_Scenarios=Ensemble_Scenarios, 
                                            imts=imts, msr=msr, rupture_aratio=rupture_aratio, 
                                            tectonicRegionType=tectonicRegionType, context_maker=context_maker, 
                                            Site_Collection=Site_Collection, correlation_model=correlation_model, 
                                            crosscorr_model=crosscorr_model, gmpes=gmpes, Weighted_Num_Realiz=Weighted_Num_Realiz)
                    chunk_results.append(result)

                results.extend(chunk_results)

            pool.close()
            pool.join()    

        # Combine results
        for result in results:
            k = result[0]
            gmf = result[1]
            for m in range(len(imts)):
                GMPEsRealizationsForProbShakeMap_AllGMPEs[j][k][m] = gmf[m]

                # PRINTING INFO
                for g, gmpe in enumerate(gmpes):    
                    if j == 0 and k == 0:
                        # Print this for one source scenario only
                        print("IMT: ", imts[m], "-- GMPE", gmpe, "is sampled", Weighted_Num_Realiz[g], "times over a total of", NumGMPEsRealizations, "times")
                    
    print("********* AGGREGATE RESULTS *******")
    # Aggregate the generated gmf at each site for Probabilistic Shakemap

    # Structure of GMPEsRealizationsForProbShakeMap 

    # 1st Index: Ensemble index (one of the N randomly sampled ensemble)
    # 2nd Index: Scenario index
    # 3rd Index: Intensity Measure index 
    # 4th Index: GMPE index
    # 5th Index: Site index 

    # For each site, there are as many values as Weighted_Num_Realiz[m], i.e. the number of realizations for the current GMPE 

    # PREPARE KEYS FOR SCENARIOS AND SITES
    keys_sites = [] 
    for s in range(len(sites)):
        keys_sites.append(f"Site_LAT:{Point(float(POIs[s][1]), float(POIs[s][0])).latitude}_LON:{Point(float(POIs[s][1]), float(POIs[s][0])).longitude}")

    keys_scen = []
    for k in range(EnsembleSize):
        keys_scen.append(f"Scenario_{k+1}") 

    # # Structure of SiteGmf 

    # # 1st Index: Index: IMT 
    # # 2nd Index: Ensemble index (one of the N randomly sampled ensemble)
    # # 3rd Index: Scenario index 
    # # 4th Index: Site index 

    # # For each site, there are as many values as NumGMPEsRealizations, i.e. tot number of realizations after merging all GMPEs
    
    #######################################
    ############# GO PARALLEL #############
    #######################################
    
    # Preallocate SiteGmf
    SiteGmf = [[[[[[] for _ in range(NumGMPEsRealizations)] for _ in range(len(sites))] for _ in range(EnsembleSize)] for _ in range(NumSampledEnsembles)] for _ in range(len(imts))]

    num_processes = num_processes 
    chunk_size_default = int(EnsembleSize/num_processes) # size of each chunk of scenarios
    print("Chunk size = ", chunk_size_default)
    last_chunk_size = chunk_size_default + EnsembleSize - num_processes * chunk_size_default # size of the last chunk
    print("Last chunk size = ", last_chunk_size)

    # Create a lock to synchronize access to the SiteGmf array
    lock = Lock()
    
    # merge the HDF5 files
    params = {} 
    for m in range(len(imts)):
        params[imts[m]] = {}
        for j in range(NumSampledEnsembles):

            params[imts[m]][j] = {}
            args = []
            for i in range(num_processes):
                if i == num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default
                args.append((m, j, i, chunk_size, num_processes, EnsembleSize, subdir1, imts, keys_scen, sites, gmpes, GMPEsRealizationsForProbShakeMap_AllGMPEs, SiteGmf, keys_sites, lock))   

            print(f"Spawning {num_processes} processes")
            processes = [Process(target=write_hdf5_file, args=arg) for arg in args]

            # start the processes
            for process in processes:
                process.start()
            for process in processes:
                process.join()

            print("Merging chunks")
            chunk_filenames = []
            merged_filename = subdir1 + "/" + "SIZE_%d" % EnsembleSize +  "_ENSEMBLE_%d_%s.hdf5" % (j + 1, imts[m])
            print("Final output .HDF5 file = ", merged_filename)
            with h5py.File(merged_filename, 'w') as merged_file:
                for i in range(num_processes):
                    chunk_filename = subdir1 + "/" + "SIZE_%d" % EnsembleSize + "_ENSEMBLE_%d_%s_chunk%d-%d.hdf5" % (j + 1, imts[m], i, i+1)
                    with h5py.File(chunk_filename, 'r') as chunk_file:
                        for scen_key in chunk_file.keys():
                            merged_file.copy(chunk_file[scen_key], scen_key)
                    chunk_filenames.append(chunk_filename) 

            # Remove all of the chunk files at once
            for chunk_filename in chunk_filenames:
                os.remove(chunk_filename)     

            # Fill params dict
            params[imts[m]][j]['Lat_Event'] = Lat_Event
            params[imts[m]][j]['Lon_Event'] = Lon_Event
            params[imts[m]][j]['event_dir'] = event_dir
            params[imts[m]][j]['NumGMPEsRealizations'] = NumGMPEsRealizations 
            params[imts[m]][j]['HDF5_Filename'] = "SIZE_%d" % EnsembleSize +  "_ENSEMBLE_%d_%s.hdf5" % (j + 1, imts[m])            

    # # CHECK
    # print(len(SiteGmf))
    # print(len(SiteGmf[0]))
    # print(len(SiteGmf[0][0]))
    # print(len(SiteGmf[0][0][0]))
    # print(len(SiteGmf[0][0][0][0]))

    print("********* DONE! *******")
    print("--- %s seconds ---" % (time.time() - start_time))

    return params

if __name__ == '__main__':
    run_prob_analysis()

