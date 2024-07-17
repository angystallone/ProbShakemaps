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

import shutil
import numpy as np
from numpy import save
import os
import importlib

import xml.etree.ElementTree as ET
import h5py

from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colorbar import cm

import multiprocessing
from multiprocessing import Pool, Process, Lock, Manager
import itertools

from functools import partial
import config

import random
from scipy import constants
import json
import sys
import matplotlib.patches as mpatches
import seaborn as sns
from openquake.hazardlib.geo import Point
from pyproj import Geod
geod = Geod(ellps='WGS84')

# FUNCTION UTILITIES

def dist_lonlat(lon1,lat1,lon2,lat2,coordtype):

    """
    Find dist (km) btw 2 points given their lon,lat
    """
    
    R = 6371.0
    
    if coordtype == 'degree':
        lon1 = lon1 * constants.pi / 180
        lat1 = lat1 * constants.pi / 180
        lon2 = lon2 * constants.pi / 180
        lat2 = lat2 * constants.pi / 180
        
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2))**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    dist_lonlat = R * c
    
    return dist_lonlat

    
def weighted_percentile(data, weights, perc):

    """
    perc : percentile in [0-1]!
    """
    ix = np.argsort(data)
    data = data[ix] # sort data
    weights = weights[ix] # sort weights
    cdf = (np.cumsum(weights) - 0.5 * weights) / np.sum(weights) # 'like' a CDF function
    return np.interp(perc, cdf, data)


def get_pois_coordinates_from_file(path, POIs_File):

    """
    Extract POIs coords from POIs file 
    """

    file = path + POIs_File

    POIs_coordinates = []
    with open(file, 'r') as f:
        for line in f:
            lat = float(line.strip().split()[0])
            lon = float(line.strip().split()[1])
            POIs_coordinates.append((lat, lon))

    return POIs_coordinates  


def get_pois(POIs_File):

    """
    Returns POIs coords and names + number of POIs in the file
    """

    path = os.path.join(os.getcwd(), "INPUT_FILES/")
    POIs_coordinates = get_pois_coordinates_from_file(path, POIs_File)
    LATs = [tup[0] for tup in POIs_coordinates]
    LONs = [tup[1] for tup in POIs_coordinates]

    POIs_lat, POIs_lon, POIs_NAMES = [], [], []

    for lat, lon in zip(LATs, LONs):

        POIs_lat.append(lat)
        POIs_lon.append(lon)
        POIs_NAMES.append(f"Site_LAT:{Point(float(lon),float(lat)).latitude}_LON:{Point(float(lon),float(lat)).longitude}")
        n_pois = len(POIs_NAMES)

    return POIs_lat, POIs_lon, POIs_NAMES, n_pois     
        

def get_pois_subset(POIs_File, Lon_Event, Lat_Event, pois_selection_method, n_pois, max_distance):

    """
    Extracts POIs subset. Two options available: 
    1. 'random': randomly extracts the subset
    2. 'azimuth_uniform': extracts POIs that are azimuthally uniformly distributed
    """

    path = os.path.join(os.getcwd(), "INPUT_FILES/")
    POIs_coordinates = get_pois_coordinates_from_file(path, POIs_File)
    LATs = [tup[0] for tup in POIs_coordinates]
    LONs = [tup[1] for tup in POIs_coordinates]

    print(f"Max distance for POIs selection = {max_distance} km")
    print("POIs selection method = ", pois_selection_method)

    if pois_selection_method == 'random':

        random_indices = random.sample(range(len(LONs)), len(LONs))
        random_lon = [LONs[i] for i in random_indices]
        random_lat = [LATs[i] for i in random_indices]

        idx_POIs, POIs_lat, POIs_lon, POIs_NAMES = [], [], [], []

        count = 0
        for idx, lon, lat in zip(random_indices, random_lon, random_lat):

            dist = dist_lonlat(lon, lat, Lon_Event, Lat_Event, 'degree')

            if dist < max_distance:

                idx_POIs.append(idx)
                POIs_lat.append(lat)
                POIs_lon.append(lon)
                POIs_NAMES.append(f"Site_LAT:{Point(float(lon),float(lat)).latitude}_LON:{Point(float(lon),float(lat)).longitude}")      
                count += 1
                if count >= n_pois:
                    break        
                    
        if count < n_pois:
            print(f"Number of POIs found = {count}")
            print("Consider changing max_distance value")
            sys.exit()       
        else:
            print("Extracted ", count, "POIs")                

        path = os.path.join(os.getcwd(), "OUTPUT")
        pois_file = os.path.join(path, "POIs.txt")
        with open(pois_file, "w") as f:
            for idx, lat, lon in zip(idx_POIs, POIs_lat, POIs_lon):
                poi = lat, lon
                f.write("{} {:.6f} {:.6f}\n".format(idx, *poi))             

        return idx_POIs, POIs_lat, POIs_lon, POIs_NAMES   

    if pois_selection_method == 'azimuth_uniform':

        DIST_BUFFER = max_distance/10
        indices = list(range(len(LONs)))

        azimuths = []
        for i in range(len(LONs)):
            azimuth, _, _ = geod.inv(LATs[i], LONs[i], Lat_Event, Lon_Event)
            azimuths.append((azimuth + 180) % 360)  # wrap azimuths to [0, 360)

        # Divide the azimuth range into equal segments and select the POI closest to the center of each segment
        segment_size = 360 / 4
        segment_centers = np.arange(segment_size/2, 360, segment_size)

        idx_POIs, POIs_lat, POIs_lon, POIs_NAMES, azimuth_POIs = [], [], [], [], []
        found_poi_count = 0
        pois_per_segment = int(n_pois/4)
        print("pois_per_segment =", pois_per_segment)
        for center in segment_centers:
            segment_indices = [j for j in indices if center - segment_size/2 <= azimuths[j] < center + segment_size/2]
            random.shuffle(segment_indices) 
            segment_pois = []
            for j in segment_indices:
                azimuth_diff = abs(azimuths[j] - center)
                if azimuth_diff < segment_size:
                    distance = dist_lonlat(LONs[j], LATs[j], Lon_Event, Lat_Event, 'degree')
                    if max_distance - DIST_BUFFER <= distance <= max_distance + DIST_BUFFER:
                        segment_pois.append(j)
                        idx_POIs.append(j)
                        POIs_lat.append(LATs[j])
                        POIs_lon.append(LONs[j])
                        POIs_NAMES.append(f"Site_LAT:{Point(float(LONs[j]), float(LATs[j])).latitude}_LON:{Point(float(LONs[j]), float(LATs[j])).longitude}")
                        azimuth_POIs.append(azimuths[j])
                        if len(segment_pois) == pois_per_segment:
                            break

            if len(segment_pois) < pois_per_segment:
                found_poi_count += 1
                print(f"WARNING! No POIs have been found between {center - segment_size/2}° and {center + segment_size/2}°") 
            
        if found_poi_count != 0:
            print(f"Number of POIs found = {found_poi_count}")
            print("Consider changing max_distance value")
            sys.exit()       
        else:
            print("Extracted ", len(azimuth_POIs), "POIs")

        # Sort idxs, lons and lats by ascending order of the azimuths
        comb_lists = list(zip(idx_POIs, azimuth_POIs))
        sorted_lists = sorted(comb_lists, key=lambda x: x[1])
        sorted_idx_POIs = [item[0] for item in sorted_lists]
        sorted_POIs_lon = [LONs[i] for i in sorted_idx_POIs]
        sorted_POIs_lat = [LATs[i] for i in sorted_idx_POIs]
        sorted_azimuths = [item[1] for item in sorted_lists]
        sorted_azimuths = [int(az) for az in sorted_azimuths]

        path = os.path.join(os.getcwd(), "OUTPUT")
        pois_file = os.path.join(path, "POIs.txt")
        with open(pois_file, "w") as f:
            for idx, lat, lon, azim in zip(sorted_idx_POIs, sorted_POIs_lat, sorted_POIs_lon, sorted_azimuths):
                poi = lat, lon
                f.write("{} {:.6f} {:.6f} {}\n".format(idx, *poi, azim))  

        return sorted_idx_POIs, sorted_POIs_lat, sorted_POIs_lon, POIs_NAMES, sorted_azimuths
    

def share_pois(POIs_File):

    """
    Shares the same POIs subset across the prob_tools
    """

    start_path = os.path.join(os.getcwd(), "OUTPUT/")
    end_path = os.path.join(os.getcwd(), "INPUT_FILES/")
    shutil.copy2(os.path.join(os.getcwd(), start_path, POIs_File), os.path.join(os.getcwd(), end_path, POIs_File))

    file_path = os.path.join(end_path, POIs_File)
    
    # Determine the method based on the first line of the file
    with open(file_path, 'r') as f:
        line = f.readline().strip()
        line_parts = line.split()
        method = 'random' if len(line_parts) == 3 else 'azimuth_uniform'
    
    idx_POIs, POIs_lat, POIs_lon, POIs_NAMES = [], [], [], []
    azimuths = None
    
    with open(file_path, 'r') as f:
        if method == 'random':
            for line in f:
                idx, lat, lon = line.strip().split()
                idx_POIs.append(int(idx))
                POIs_lat.append(float(lat))
                POIs_lon.append(float(lon))
                POIs_NAMES.append(f"Site_LAT:{Point(float(lon), float(lat)).latitude}_LON:{Point(float(lon), float(lat)).longitude}")
            n_pois = len(POIs_NAMES)
            return idx_POIs, POIs_lat, POIs_lon, POIs_NAMES, azimuths, n_pois
        else:
            azimuths = []
            for line in f:
                idx, lat, lon, azimuth = line.strip().split()
                idx_POIs.append(int(idx))
                POIs_lat.append(float(lat))
                POIs_lon.append(float(lon))
                azimuths.append(float(azimuth))
                POIs_NAMES.append(f"Site_LAT:{Point(float(lon), float(lat)).latitude}_LON:{Point(float(lon), float(lat)).longitude}")
            n_pois = len(POIs_NAMES)
            return idx_POIs, POIs_lat, POIs_lon, POIs_NAMES, azimuths, n_pois    


def pois_map(POIs_lat, POIs_lon, Lat_Event, Lon_Event, deg_round, path):

    """
    Returns a map with POIs and event
    """

    poi_indices = [idx + 1 for idx in range(len(POIs_lat))]

    xlim_min = np.min(np.floor(POIs_lon/deg_round) * deg_round)
    xlim_max = np.max(np.ceil(POIs_lon/deg_round) * deg_round)
    ylim_min = np.min(np.floor(POIs_lat/deg_round) * deg_round)
    ylim_max = np.max(np.ceil(POIs_lat/deg_round) * deg_round)


    latitudes = np.arange(-90, 91, 2)
    longitudes = np.arange(-180, 181, 2)

    m = Basemap(projection='merc',llcrnrlat=ylim_min,urcrnrlat=ylim_max,\
                llcrnrlon=xlim_min,urcrnrlon=xlim_max,lat_ts=20,resolution='i') 

    m.drawcoastlines(linewidth=0.5)   
    m.drawstates(linewidth=0.5)

    m.drawparallels(latitudes, labels=[1,0,0,0], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])
    m.drawmeridians(longitudes, labels=[0,0,0,1], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])

    m.drawmapboundary(linewidth=2, color='black', fill_color='white')

    x, y = m(POIs_lon, POIs_lat)
    x_event, y_event = m(Lon_Event, Lat_Event)

    for i, label in enumerate(poi_indices):
        plt.text(x[i], y[i] + 8000, label, fontsize=5, color='black')
    
    m.scatter(x, y, s=10, marker='o', color='red', label="POIs")

    m.scatter(x_event, y_event, s=70, marker='*', color='blue', label="Epicenter")
    plt.legend(loc='lower left')

    plt.savefig(path + "/POIs_Map.pdf", dpi=300, bbox_inches='tight')
    plt.close()


def get_params():
        
    """
    Returns params needed for prob analysis
    """

    config_dict = config.load_config('input_file.txt')
    ID_Event = config_dict['ID_Event']

    # Install dir and event dir
    install_path, data_path = get_config_paths()
    event_dir = os.path.join(data_path, ID_Event, "current")

    if not os.path.exists(event_dir):
        raise NotADirectoryError(f"{event_dir} is not a valid directory.")

    # Load and parse the event file
    eventxml = os.path.join(event_dir, "event.xml")    
    if not os.path.isfile(eventxml):
        raise FileNotFoundError(f"{eventxml} does not exist.")   

    tree = ET.parse(eventxml)
    root = tree.getroot()

    # Get the latitude and longitude of the event
    Lat_Event = root.attrib.get('lat') 
    Lon_Event = root.attrib.get('lon') 

    listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
    scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]

    # Get the number of scenarios in the file
    with open(os.path.join(listscenarios_dir, scenarios_file[0]), 'r') as f:
        EnsembleSize = 0
        for _ in f:
            EnsembleSize += 1

    params = {}
    params['install_path'] = install_path
    params['data_path'] = data_path
    params['Lat_Event'] = Lat_Event
    params['Lon_Event'] = Lon_Event
    params['event_dir'] = event_dir
    params['Ensemble_Size'] = EnsembleSize

    return params
    

##############################################################################
############################ PROBSHAKEMAP ####################################
##############################################################################

class Main:
    def __init__(self, IMT, pois_file, NumGMPEsRealizations, num_processes):

        self.imt = IMT
        self.pois_file = pois_file
        self.NumGMPEsRealizations = NumGMPEsRealizations
        self.num_processes = num_processes

        params = get_params()
        self.event_dir = params["event_dir"]
        self.EnsembleSize = params['Ensemble_Size']
        self.install_path = params['install_path'] 
        self.data_path = params['data_path'] 
        
        self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
        self.POIs_lat = np.array(self.POIs_lat)
        self.POIs_lon = np.array(self.POIs_lon)

    def process_scenario(scen, Ensemble_Scenarios, msr, rupture_aratio, 
                        tectonicRegionType, context_maker, Site_Collection, 
                        correlation_model, crosscorr_model, gmpes, Weighted_Num_Realiz):
        
        """
        For a given scenario, retrieves GMFs from all GMPEs at all POIs
        """

        # Get scenario index 
        k = Ensemble_Scenarios.index(scen)

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

        # Get an array of shape (4, G, M, N) with mean and stddevs 
        # (G = number of GMPEs, M = number of IMTs, N = number of sites)
        mean_and_stdev = context_maker.get_mean_stds(ctx)  
            
        # Loop over GMPEs    
        gmf = []
        for g, gmpe in enumerate(gmpes):
            # 'gc.compute' --> Compute gmf and returns an array of shape (num_imts, num_sites, num_events) with the sampled gmf, and two arrays with shape 
            # (num_imts, num_events): sig for tau and eps for the random part
            # Note that Weighted_Num_Realiz[g] == num_events
            gf = gc.compute(gmpe, Weighted_Num_Realiz[g], mean_and_stdev[:, g, :, :])  

            # Append only first array output from gc.compute (shape: (num_imts, num_sites, num_events))  
            gmf.append(gf[0])

        return k, gmf
    
    def aggregate_gmfs(scen, Ensemble_Scenarios, NumGMPEsRealizations, sites, gmpes_list, GMPEsRealizationsForProbShakeMap_AllGMPEs):
        
        """
        For a given scenario, aggregates GMFs from all GMPEs at each POI
        """

        # Get scenario index 
        k = Ensemble_Scenarios.index(scen)

        SiteGmf_scen = np.empty((len(sites), NumGMPEsRealizations), dtype=object)

        # Loop over sites
        for s in range(len(sites)):
            SiteGmfGMPE = []
            for g in range(len(gmpes_list)): 
                # IMT index = 0 as we consider only 1 IMT at a time
                SiteGmfGMPE.append(GMPEsRealizationsForProbShakeMap_AllGMPEs[k][g][0][s])

            SiteGmf_scen[s] = [x for sublist in SiteGmfGMPE for x in sublist]

        return k, SiteGmf_scen
    
    def run_prob_analysis(self):

        """
        Runs the prob analysis
        """

        print("********* STARTING PROB ANALYSIS *******")

        # Load configuration parameters
        config_dict = config.load_config('input_file.txt')

        tectonicRegionType = config_dict['tectonicRegionType']
        mag_scaling = config_dict['mag_scaling']
        rupture_aratio = config_dict['rupture_aratio']
        ID_Event = config_dict['ID_Event']
        vs30file = config_dict['vs30file']
        CorrelationModel = config_dict['CorrelationModel']
        CrosscorrModel = config_dict['CrosscorrModel']
        vs30_clustering = config_dict['vs30_clustering']
        truncation_level = config_dict['truncation_level']
        seed = config_dict['seed']

        print("Install Path = ", self.install_path)
        print("Data Path = ", self.data_path)
        print("Number of source scenarios to process = ", self.EnsembleSize)
        print("Number of CPU processes: ", str(self.num_processes))

        path = os.path.join(os.getcwd(), "OUTPUT")
        if not os.path.exists(path):
            os.makedirs(path)

        # PRINT USER'S INPUT 
        print("TectonicRegionType: " + tectonicRegionType)
        print("Importing " + mag_scaling + " as magnitude scaling relationship")
        module = importlib.import_module('openquake.hazardlib.scalerel')
        msr = getattr(module, mag_scaling)
        print("Rupture aspect ratio: " + str(rupture_aratio))
        print("Event ID: " + ID_Event)
        print("POIs file: " + self.pois_file)
        if vs30file == None:
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
        print("Seed: " + str(seed))
        print("Number of GMPEs realizations per POI: " + str(self.NumGMPEsRealizations))
        print("Intensity measure type: " + str(self.imt))

        # Collects event and configuration data and creates the file shake_data.hdf
        assemble = AssembleModule(ID_Event, comment='Test comment.')
        assemble.execute()

        # Reads the data in shake_data.hdf and produces an interpolated ShakeMap --> shake_result.hdf
        model = ModelModule(ID_Event)
        model.execute()

        # Generate model_select.conf file, containing GMPE sets 
        select = SelectModule(ID_Event)
        select.execute()

        print("********* RETRIEVING GMPEs *******")
        
        # Read the event.xml file and generate a GMPE set for the event based on the event’s residence within, 
        # and proximity to, a set of predefined tectonic regions and user-defined geographic areas

        # Extract GMPE set
        conf_filename = self.event_dir + '/model_select.conf'
        config_model_select = ConfigObj(conf_filename)
        print("Config filename = ", conf_filename)

        GMPE_Set = []
        for key, _ in config_model_select['gmpe_sets'].items():
            GMPE_Set.append(config_model_select['gmpe_sets'][key]['gmpes'][0])

        print("GMPEs Set available = ", GMPE_Set)

        # config files
        conf_filename = self.install_path + '/config/gmpe_sets.conf'
        config_gmpe_sets = ConfigObj(conf_filename)

        conf_filename = self.install_path + '/config//modules.conf'
        config_modules = ConfigObj(conf_filename)

        # Get GMPEs acronyms
        GMPEs_Weights = {}
        for key, _ in config_gmpe_sets['gmpe_sets'].items():
            if key in GMPE_Set:
                if not isinstance(config_gmpe_sets['gmpe_sets'][key]['gmpes'], list):
                    acronym = str(config_gmpe_sets['gmpe_sets'][key]['gmpes'])
                    weight = float(config_gmpe_sets['gmpe_sets'][key]['weights'])
                    GMPEs_Weights[acronym] = weight
                else:
                    for acronym, weight in zip(config_gmpe_sets['gmpe_sets'][key]['gmpes'],
                                               config_gmpe_sets['gmpe_sets'][key]['weights']):
                        GMPEs_Weights[str(acronym)] = float(weight)

        # Get GMPEs from the acronyms
        GMPEs_Names = {}
        for key, item in config_modules['gmpe_modules'].items():  
            if key in GMPEs_Weights.keys():
                print(f"Importing {item[0]}")
                GMPEs_Names[key] = item[0]

        # Convert gmpes names into GSIM instances
        gmpes = {}
        for acronym, elem in GMPEs_Names.items():
            gmpes[acronym] = valid.gsim(elem)  # OpenQuake equivalent of getattr

        # Filter GMPEs with the selected IMT available
        gmpes_ok = {}
        for acronym, gmpe in gmpes.items():
            list_of_imts = ', '.join([imt.__name__ for imt in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES])   
            if self.imt in ['PGA', 'PGV']:
                if self.imt in list_of_imts:
                    gmpes_ok[acronym] = gmpe
            else:
                # SA
                if self.imt[:2] in list_of_imts:
                    gmpes_ok[acronym] = gmpe

        gmpes = gmpes_ok

        if len(gmpes) != 0:
            print("GMPEs with the requested IMT available = ", [item for _, item in gmpes.items()])     
        else:
            print('No GMPEs with the requested IMT available')   
            sys.exit() 

        # Update GMPEs_Names and GMPEs_Weights
        GMPEs_Weights = {acronym: weight for acronym, weight in GMPEs_Weights.items() if acronym in gmpes.keys()}
        GMPEs_Names = {acronym: name for acronym, name in GMPEs_Names.items() if acronym in gmpes.keys()}   

        if vs30file is not None:
            print("********* LOADING Vs30 *******")
            vs30fullname = os.path.join(self.data_path, 'shakemap_data', 'vs30', vs30file)
            vs30grid = GMTGrid.load(vs30fullname)
            # Interpolate Vs30 values at POIs 
            vs30_POIs = vs30grid.getValue(self.POIs_lat, self.POIs_lon, method="nearest")

        print("********* DEFINING OpenQuake SITE COLLECTION *******")

        # Define a SiteCollection for all the POIs
        sites = []
        for i in range(len(self.POIs_NAMES)):
            site_location = Point(self.POIs_lon[i], self.POIs_lat[i])
            if vs30file == None:
                # If Vs30 file is not provided, use default value for Vs30
                site = Site(location=site_location, vs30=760., vs30measured=False, z1pt0=40., z2pt5=1.0)
            else:    
                site = Site(location=site_location, vs30=vs30_POIs[i], vs30measured=False, z1pt0=40., z2pt5=1.0)
            sites.append(site)
            
        Site_Collection = SiteCollection(sites)
        print(Site_Collection.complete)    

        print("********* BUILDING OPENQUAKE CONTEXTS *******")

        # Build OpenQuake contexts

        # # Define input parameters for ContextMaker

        imtls = {}
        imtls[self.imt] = []
        param = dict(imtls=imtls)   

        # Instantiate a ContextMaker object (Note: independent from source and sites!)
        gmpes_list = list(gmpes.values())
        context_maker = ContextMaker(tectonicRegionType, gmpes_list, param)

        print("********* SAMPLING UNCERTAINTY *******")

        correlation_model = correlation_model(vs30_clustering=vs30_clustering)
        crosscorr_model = crosscorr_model(truncation_level=truncation_level)

        # This is needed to sample proportionally to the weight of the GMPEs
        Weighted_Num_Realiz = []
        for acronym, weight in GMPEs_Weights.items():
            Weighted_Num_Realiz.append(round(self.NumGMPEsRealizations * GMPEs_Weights[acronym])) 

        if 0 in Weighted_Num_Realiz:
            raise RuntimeError("Increase NumGMPEsRealizations to sample all the GMPEs")   

        current_total_samples = sum(Weighted_Num_Realiz)
        
        if current_total_samples < self.NumGMPEsRealizations:
            remaining_samples = self.NumGMPEsRealizations - current_total_samples
            while remaining_samples > 0:
                # Assign remaining samples to the GMPE with the highest weight 
                GMPEs_Weights_list = list(GMPEs_Weights.values())
                max_weight_index = GMPEs_Weights_list.index(max(GMPEs_Weights_list))
                Weighted_Num_Realiz[max_weight_index] += 1
                remaining_samples -= 1

        total_assigned = sum(Weighted_Num_Realiz)
        # print(total_assigned)
 
        # Sample from the total variability of ground motion taking into account both inter- and intra-event variability (for one source scenario only)
        # gmf = exp(mu + crosscorel(tau) + spatialcorrel(phi)) --> See: https://docs.openquake.org/oq-engine/advanced/latest/event_based.html#correlation-of-ground-motion-fields

        GMPEsRealizationsForProbShakeMap_AllGMPEs = [None] * self.EnsembleSize

        listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
        scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]
        f = open(os.path.join(listscenarios_dir, scenarios_file[0]), 'r')
        print("List of scenarios: ", scenarios_file[0])

        Ensemble_Scenarios = []
        for k, line in enumerate(f):
            scen = line.strip().split(' ')
            Ensemble_Scenarios.append(scen)  

        # Set the random seed for reproducibility in OpenQuake GmfComputer
        np.random.seed(seed)

        ################################
        # SETTING MULTIPROCESSING PARAMS
        ################################

        chunk_size_default = int(self.EnsembleSize/self.num_processes) # size of each chunk of scenarios
        #print("Chunk size = ", chunk_size_default)
        last_chunk_size = chunk_size_default + self.EnsembleSize - self.num_processes * chunk_size_default # size of the last chunk
        #print("Last chunk size = ", last_chunk_size)

        # Create pool of worker processes
        with Pool(processes=self.num_processes) as pool:
            results = []

            # iterate over processes
            for i in range(self.num_processes):
                if i == self.num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default

                start_idx = i * chunk_size
                end_idx = (i+1) * chunk_size
                # adjust k_start and k_end for the last chunk
                if i == self.num_processes - 1:
                    start_idx = self.EnsembleSize - chunk_size
                    end_idx = self.EnsembleSize 

                chunk = Ensemble_Scenarios[start_idx:end_idx]

                chunk_results = []
                for scenario in chunk:
                    result = Main.process_scenario(scen=scenario, Ensemble_Scenarios=Ensemble_Scenarios, 
                                            msr=msr, rupture_aratio=rupture_aratio, tectonicRegionType=tectonicRegionType,
                                            context_maker=context_maker, Site_Collection=Site_Collection, 
                                            correlation_model=correlation_model, crosscorr_model=crosscorr_model, gmpes=gmpes_list, 
                                            Weighted_Num_Realiz=Weighted_Num_Realiz)
                    chunk_results.append(result)

                results.extend(chunk_results)

            pool.close()
            pool.join()    

        # Combine results
        for result in results:
            k = result[0]
            gmf = result[1]
            GMPEsRealizationsForProbShakeMap_AllGMPEs[k] = gmf

            # PRINTING INFO
            for g, gmpe in enumerate(gmpes_list):    
                if k == 0:
                    # Print this for one source scenario only
                    print("IMT: ", self.imt, "-- GMPE", gmpe, "is sampled", Weighted_Num_Realiz[g], "times over a total of", self.NumGMPEsRealizations, "times")

        # GMFs AGGREGATION
        # Aggregate the generated gmf at each site for Probabilistic Shakemap

        # Structure of GMPEsRealizationsForProbShakeMap_AllGMPEs
        # 1st Index: Scenario index
        # 2nd Index: GMPE index
        # 3rd Index: IMT index
        # 4th Index: Site index 

        # CHECK
        # print("SHAPE = ", len(GMPEsRealizationsForProbShakeMap_AllGMPEs))
        # print("SHAPE = ", len(GMPEsRealizationsForProbShakeMap_AllGMPEs[0]))
        # print("SHAPE = ", len(GMPEsRealizationsForProbShakeMap_AllGMPEs[0][0]))
        # print("SHAPE = ", len(GMPEsRealizationsForProbShakeMap_AllGMPEs[0][0][0]))

        # For each site, there are as many values as the number of realizations for the current GMPE 

        # PREPARE KEYS FOR SCENARIOS AND SITES
        keys_sites = [] 
        for s in range(len(sites)):
            keys_sites.append(f"Site_LAT:{Point(float(self.POIs_lon[s]), float(self.POIs_lat[s])).latitude}_LON:{Point(float(self.POIs_lon[s]), float(self.POIs_lat[s])).longitude}")

        keys_scen = []
        for k in range(self.EnsembleSize):
            keys_scen.append(f"Scenario_{k+1}") 

        # # Structure of SiteGmf 
        # 1st Index: Scenario index 
        # 2nd Index: Site index 

        # Preallocate SiteGmf
        SiteGmf = np.empty((self.EnsembleSize, len(sites), self.NumGMPEsRealizations), dtype=object)

        # Create pool of worker processes
        with Pool(processes=self.num_processes) as pool:
            results = []

            # iterate over processes
            for i in range(self.num_processes):
                if i == self.num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default

                start_idx = i * chunk_size
                end_idx = (i+1) * chunk_size
                # adjust k_start and k_end for the last chunk
                if i == self.num_processes - 1:
                    start_idx = self.EnsembleSize - chunk_size
                    end_idx = self.EnsembleSize 

                chunk = Ensemble_Scenarios[start_idx:end_idx]

                chunk_results = []
                for scenario in chunk:
                    result = Main.aggregate_gmfs(scen=scenario, Ensemble_Scenarios=Ensemble_Scenarios, 
                                                 NumGMPEsRealizations=self.NumGMPEsRealizations, sites=sites, gmpes_list=gmpes_list, 
                                                 GMPEsRealizationsForProbShakeMap_AllGMPEs=GMPEsRealizationsForProbShakeMap_AllGMPEs)
                    chunk_results.append(result)

                results.extend(chunk_results)

            pool.close()
            pool.join()    

        # Combine results
        for result in results:
            i_scen = result[0]
            SiteGmf_scen = result[1]
            SiteGmf[i_scen] = SiteGmf_scen

        print("********* PROB ANALYSIS DONE! *******")

        # CHECK
        # print(SiteGmf_scen.shape)
        # print(SiteGmf.shape)

        prob_output = {
            "SiteGmf": SiteGmf,
            "keys_scen": keys_scen,
            "keys_sites": keys_sites
        }

        return prob_output
    

class Write():
    def __init__(self, IMT, EnsembleSize, keys_scen, SiteGmf, keys_sites, num_processes):

        self.num_processes = num_processes
        self.imt = IMT
        self.EnsembleSize = EnsembleSize
        self.keys_scen = keys_scen
        self.SiteGmf = SiteGmf
        self.keys_sites = keys_sites

    def write_output(self):

        print("Numer of POIs = ", len(self.keys_sites), "POIs")

        path = os.path.join(os.getcwd(), "OUTPUT")
        if not os.path.exists(path):
            os.makedirs(path)

        subdir1 = path + "/HDF5_FILES" 
        if not os.path.exists(subdir1):
            os.makedirs(subdir1)

        filename = os.path.join(subdir1, "SIZE_%d" % self.EnsembleSize +  "_ENSEMBLE_%s.hdf5" % (self.imt))

        print("********* WRITING FILE *******")

        # One .hdf5 file per ensemble
        h5f = h5py.File(filename,'w')

        for k in range(self.EnsembleSize):

            # One group per scenario
            grp = h5f.create_group(self.keys_scen[k])

            # One dataset per site
            for s in range(len(self.keys_sites)):    
        
                # WRITE OUTPUT FILE    
                data_name = self.keys_sites[s]
                grp.create_dataset(data_name, data=self.SiteGmf[k][s], compression="gzip", compression_opts=9)
        
        h5f.close()

        print("********* OUTPUT FILE READY! *******")


class StationRecords:
    def __init__(self, IMT, imt_min, imt_max, deg_round, stationfile):
                 
        self.imt = IMT
        self.imt_min = imt_min
        self.imt_max = imt_max
        self.deg_round = deg_round 
        self.stationfile = stationfile

        params = get_params()
        self.event_dir = params["event_dir"]

    def get_data_coord(self):    

        file_station = os.path.join(self.event_dir, self.stationfile)

        with open(file_station) as json_file:

            data = json.load(json_file)
            data_lon = np.zeros((len(data['features'][:])))
            data_lat = np.zeros((len(data['features'][:])))
            for i in range(len(data['features'][:])):
                data_lon[i] = data['features'][i]['geometry']['coordinates'][0]
                data_lat[i] = data['features'][i]['geometry']['coordinates'][1]

        return data_lon, data_lat
    
    def get_data(self):    

        file_station = os.path.join(self.event_dir, self.stationfile)

        with open(file_station) as json_file:

            key = self.imt.lower()
            data = json.load(json_file)
            data_imt = np.zeros((len(data['features'][:])))
            for i in range(len(data['features'][:])):
                if self.imt in ['PGA', 'PGV']:
                    tmp = data['features'][i]['properties'][key]
                else:
                    list_predictions = data['features'][i]['properties']['predictions']
                    for pred in list_predictions:
                        if pred['name'] == key:
                            tmp = pred['value']
                if tmp != 'null':
                    if self.imt == 'PGV':
                        data_imt[i] = tmp 
                    else:
                        data_imt[i] = tmp / 100 
                else:
                    data_imt[i] = -0.01
                        
        print(f"{self.imt} range: {np.min(data_imt)} - {np.max(data_imt)}")

        return data_imt
    
    def get_stations(self):    

        file_station = os.path.join(self.event_dir, self.stationfile)

        with open(file_station) as json_file:

            data = json.load(json_file)
            station_id = [] * len(data['features'][:])
            station_name = [] * len(data['features'][:])
            for i in range(len(data['features'][:])):
                station_id.append(data['features'][i]['properties']['code'])
                station_name.append(data['features'][i]['properties']['name'])

        return station_id, station_name

    def plot(self):

        print("********* PLOTTING STATIONS DATA *******")

        data_lon, data_lat = StationRecords.get_data_coord(self)
        data_imt = StationRecords.get_data(self)
  
        xlim_min = np.min(np.floor(data_lon/self.deg_round) * self.deg_round)
        xlim_max = np.max(np.ceil(data_lon/self.deg_round) * self.deg_round)
        ylim_min = np.min(np.floor(data_lat/self.deg_round) * self.deg_round)
        ylim_max = np.max(np.ceil(data_lat/self.deg_round) * self.deg_round)

        latitudes = np.arange(-90, 91, 2)
        longitudes = np.arange(-180, 181, 2)

        cm = plt.cm.get_cmap('YlOrRd')

        fig= plt.figure(figsize=(9, 6))
        m = Basemap(projection='merc',llcrnrlat=ylim_min,urcrnrlat=ylim_max,\
                    llcrnrlon=xlim_min,urcrnrlon=xlim_max,lat_ts=20,resolution='i')
        m.drawcoastlines()

        m.drawparallels(latitudes, labels=[1,0,0,0], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])
        m.drawmeridians(longitudes, labels=[0,0,0,1], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])

        m.drawmapboundary(linewidth=2, color='black', fill_color='white')

        x, y = m(data_lon, data_lat)

        plt.title(f"{self.stationfile}, {self.imt}")
        sc = m.scatter(x, y, c=np.log10(data_imt), vmin=np.log10(self.imt_min), vmax=np.log10(self.imt_max), s=20, edgecolors='black', linewidths=0.2, cmap = cm)
        cbar = plt.colorbar(sc)
        cbar.set_label(f"log10({self.imt})")

        plt.show()

        path = os.path.join(os.getcwd(), "OUTPUT")
        if not os.path.exists(path):
            os.makedirs(path)
        figname = f"/Data_stationfile_{self.imt}.pdf"
        fig.savefig(path + figname, dpi=200)   
        print("DONE!") 
        print(f"***** Figure saved in {path} *****")


class QueryHDF5:
    def __init__(self, scenario, pois_file, pois_subset, n_pois, max_distance,
                pois_selection_method, Lon_Event, Lat_Event):

        self.scenario = scenario
        self.pois_file = pois_file
        self.pois_subset = bool(pois_subset)
        self.n_pois = n_pois
        self.max_distance = max_distance    
        self.pois_selection_method = pois_selection_method
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event

        if self.pois_subset == False:  

            self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)
            print("Found ", self.n_pois, "POIs")

        else:
            if self.pois_selection_method == 'random':
                self.azimuths = None
                self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                            self.pois_selection_method, self.n_pois, self.max_distance)
            else:
                self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                            self.pois_selection_method, self.n_pois, self.max_distance)

        outfile_dir = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/")
        filename = [name for name in os.listdir(outfile_dir) if name != ".DS_Store"][0]
        self.outputfile = os.path.join(outfile_dir, filename)
        print("File to query = ", self.outputfile)

    def get_scenarios(self):   
        with h5py.File(self.outputfile, "r") as f: 
            Scenarios = list(f)
            Scenarios_n = len(Scenarios)  

        return Scenarios, Scenarios_n

    def print_info(self):

        print("********* QUERYING FILE *******")

        _, Scenarios_n = QueryHDF5.get_scenarios(self)

        with h5py.File(self.outputfile, "r") as f: 
            scen = 'Scenario_%d' % self.scenario

            print(f"Number of Scenarios: {Scenarios_n}")
            
            data = {}
            for poi in self.POIs_NAMES:
                data[poi] = list(f[scen][poi])
            print(f"Number of GMF per POI: {len(list(f[scen][poi]))}")    
                    
            filename = os.path.join(os.getcwd(), "OUTPUT/", 'GMF_info.txt')
            with open(filename, 'w') as f:
                for poi in self.POIs_NAMES:
                    output_str = f"GMF realizations at {poi} for {scen}: {data[poi]}\n"
                    f.write(output_str)
            print(f"GMF info for {scen} printed in the file '{filename}")


class GetStatistics:
    def __init__(self, SiteGmf, EnsembleSize, Lon_Event, Lat_Event, NumGMPEsRealizations,
                  event_dir, IMT, imt_min, imt_max, fileScenariosWeights, 
                 pois_file, pois_subset, n_pois, max_distance, pois_selection_method, deg_round,
                 pois_subset_flag, num_processes, vector_npy):

        self.NumGMPEsRealizations = NumGMPEsRealizations
        self.SiteGmf = SiteGmf
        self.EnsembleSize =  EnsembleSize
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.event_dir = event_dir
        self.imt = IMT
        self.imt_min = imt_min
        self.imt_max = imt_max
        self.fileScenariosWeights = fileScenariosWeights
        self.pois_file = pois_file
        self.pois_subset = bool(pois_subset)
        self.n_pois = n_pois
        self.max_distance = max_distance
        self.pois_selection_method = pois_selection_method
        self.deg_round = deg_round
        self.pois_subset_flag = bool(pois_subset_flag)
        self.num_processes = num_processes
        self.vector_npy = bool(vector_npy)

        print(f"Save vector.npy set to {self.vector_npy}")

        if self.pois_subset_flag:
            if self.pois_subset == False:  

                self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
                self.POIs_lat = np.array(self.POIs_lat)
                self.POIs_lon = np.array(self.POIs_lon)
                print("Found ", self.n_pois, "POIs")
                self.i_pois = list(range(self.n_pois))

            else:

                if self.pois_selection_method == 'random':
                    self.azimuths = None
                    self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                                self.pois_selection_method, self.n_pois, self.max_distance)
                else:
                    self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                                self.pois_selection_method, self.n_pois, self.max_distance)
        
                self.n_pois = n_pois
                self.POIs_lat = np.array(self.POIs_lat)
                self.POIs_lon = np.array(self.POIs_lon)
                self.i_pois = self.idx_POIs

        else:
            self.pois_file = "POIs.txt"
            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths, self.n_pois = share_pois(self.pois_file)
            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)
            print("Found ", self.n_pois, "POIs in the shared file ", self.pois_file)
            self.i_pois = self.idx_POIs

    def process_scen_gmf(scen, Ensemble_Scenarios, site_gmf, i_pois, n_pois, num_realizations, weights):
        
        """
        For a given scenario, gather all GMFs at each given POI
        Needed for getting the ground motion distribution at each POI
        """

        # Get scenario index 
        i_scen = Ensemble_Scenarios.index(scen)
     
        vector_scen = np.zeros([n_pois, num_realizations])
        weight_scen = np.zeros([n_pois, num_realizations]) 

        for jp in range(n_pois):

            vector_scen[jp] = site_gmf[i_scen][i_pois[jp]]
            weight_scen = weights[i_scen] / num_realizations

        return i_scen, vector_scen, weight_scen

    def calc_statistics(self):

        """
        Calculates the statistics of the ground motion distributions at all POIs
        """

        # Define statistical measures
        vector_stat_names = ['Mean', 'Median', 'Percentile 10', 'Percentile 20', 'Percentile 80', 'Percentile 90', 
                            'Percentile 5', 'Percentile 95', 'Percentile 2_5', 'Percentile 97_5']

        # Default weights
        weights = np.ones(self.EnsembleSize) / self.EnsembleSize

        # Load weights from file if provided
        if self.fileScenariosWeights:
            print('Weights file:', self.fileScenariosWeights)
            with open(self.fileScenariosWeights) as f:
                weights = np.array([float(line.split()[0]) for line in f.readlines()])

        weights = weights / np.sum(weights)  # Normalize weights

        listscenarios_dir = os.getcwd() + "/INPUT_FILES/ENSEMBLE/"
        scenarios_file = [name for name in os.listdir(listscenarios_dir) if name != ".DS_Store"]
        f = open(os.path.join(listscenarios_dir, scenarios_file[0]), 'r')
        Ensemble_Scenarios = []
        for _, line in enumerate(f):
            scen = line.strip().split(' ')
            Ensemble_Scenarios.append(scen) 

        # Chunking
        chunk_size_default = int(self.EnsembleSize/self.num_processes) # size of each chunk of scenarios
        last_chunk_size = chunk_size_default + self.EnsembleSize - self.num_processes * chunk_size_default # size of the last chunk

        # Create pool of worker processes
        with Pool(processes=self.num_processes) as pool:
            results = []
        
            # iterate over processes
            for i in range(self.num_processes):
                if i == self.num_processes - 1:
                    chunk_size = last_chunk_size # adjust chunk size for the last process
                else:
                    chunk_size = chunk_size_default

                start_idx = i * chunk_size
                end_idx = (i+1) * chunk_size
                # adjust k_start and k_end for the last chunk
                if i == self.num_processes - 1:
                    start_idx = self.EnsembleSize - chunk_size
                    end_idx = self.EnsembleSize 

                chunk = Ensemble_Scenarios[start_idx:end_idx]

                chunk_results = []
                for scenario in chunk: 
                    result = GetStatistics.process_scen_gmf(scen=scenario, Ensemble_Scenarios=Ensemble_Scenarios, 
                                                             site_gmf=self.SiteGmf, i_pois=self.i_pois, n_pois=self.n_pois, 
                                                             num_realizations=self.NumGMPEsRealizations, weights=weights)
                    chunk_results.append(result)

                results.extend(chunk_results)

            pool.close()
            pool.join()    

        vector_stat = {name: [0] * self.n_pois for name in vector_stat_names}
        # 'vector' gathers, for each POI, the ground motion distribution collecting
        # all GMFs from all GMPEs across all scenarios 
        vector = np.zeros((self.n_pois, self.NumGMPEsRealizations * self.EnsembleSize))
        weight = np.zeros((self.n_pois, self.NumGMPEsRealizations * self.EnsembleSize)) 

        # Collect results
        for result in results:  
            i_scen = result[0]
            vector_scen = result[1]
            weight_scen = result[2]
            start = i_scen * self.NumGMPEsRealizations
            end = (i_scen + 1) * self.NumGMPEsRealizations
            vector[:, start:end] = vector_scen
            weight[:, start:end] = weight_scen

        # Calculate statistics
        for jp in range(self.n_pois):
            vector_stat['Mean'][jp] = np.sum(vector[jp, :] * weight[jp, :]) / np.sum(weight[jp, :])
            vector_stat['Median'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.5)
            vector_stat['Percentile 10'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.1)
            vector_stat['Percentile 20'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.2)
            vector_stat['Percentile 80'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.8)
            vector_stat['Percentile 90'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.9)
            vector_stat['Percentile 5'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.05)
            vector_stat['Percentile 95'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.95)
            vector_stat['Percentile 2_5'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.025)
            vector_stat['Percentile 97_5'][jp] = weighted_percentile(vector[jp, :], weight[jp, :], 0.975)

        stats = {'vector_stat': vector_stat, 'vector': vector, 'weight': weight}

        return stats, vector_stat_names

    def save_statistics(self):

        """
        Saves 'vector_stat.npy' with all statistics
        [OPTIONAL] Saves 'vector.npy' with all ground motion distributions at all POIs (can be huge!)
        """

        stats, _ = GetStatistics.calc_statistics(self)
        vector_stat = stats['vector_stat'] 
        vector = stats['vector'] 
        #weight = stats['weight'] 

        path = os.path.join(os.getcwd(), "OUTPUT/npyFiles")
        if not os.path.exists(path):
            os.makedirs(path)
        print("********* SAVING STATISTICS *******")
        save(path + '/' + f"{self.imt}" + '_vector_stat.npy', vector_stat)
        if self.vector_npy == True:
            save(path + '/' + f"{self.imt}" + '_vector.npy', vector)

    def plot_statistics(self):

        """
        Visualizes statistics from 'vector_stat.npy' on a map and saves it
        """

        print("********* GENERATING STATISTICS PLOTS *******")

        path = os.path.join(os.getcwd(), "OUTPUT/STATISTICS")
        if not os.path.exists(path):
            os.makedirs(path)

        stats, vector_stat_names = GetStatistics.calc_statistics(self)
        vector_stat = stats['vector_stat'] 

        dim_point = 10

        # TARGET POINTS
        xlim_min = np.min(np.floor(self.POIs_lon/self.deg_round) * self.deg_round)
        xlim_max = np.max(np.ceil(self.POIs_lon/self.deg_round) * self.deg_round)
        ylim_min = np.min(np.floor(self.POIs_lat/self.deg_round) * self.deg_round)
        ylim_max = np.max(np.ceil(self.POIs_lat/self.deg_round) * self.deg_round)

        for _,name in enumerate(vector_stat_names):

            if name not in ['Percentile 2_5', 'Percentile 97_5', 'Percentile 5', 'Percentile 95']:

                fig = plt.figure(figsize=(9, 6))

                latitudes = np.arange(-90, 91, 2)
                longitudes = np.arange(-180, 181, 2)

                cm = plt.cm.get_cmap('turbo')

                m = Basemap(projection='merc',llcrnrlat=ylim_min,urcrnrlat=ylim_max,\
                        llcrnrlon=xlim_min,urcrnrlon=xlim_max,lat_ts=20,resolution='i')
                m.drawcoastlines()

                m.drawparallels(latitudes, labels=[1,0,0,0], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])
                m.drawmeridians(longitudes, labels=[0,0,0,1], fontsize=8, linewidth=0.5, color='gray', dashes=[1, 2])

                m.drawmapboundary(linewidth=2, color='black', fill_color='white')

                x, y = m(self.POIs_lon, self.POIs_lat)
                
                tmp = vector_stat[name]
                sc = m.scatter(x, y, c=np.log10(tmp), vmin=np.log10(self.imt_min), vmax=np.log10(self.imt_max), s=dim_point, cmap = cm)
                plt.title(name)
                cbar = plt.colorbar(sc)
                cbar.set_label(f"log10({self.imt})")
                plt.show()

                figname = path + '/' + 'Stat_' + name.strip() + '.pdf'
                fig.savefig(figname, dpi = 200)
        print(f"***** Figures saved in {path} *****")


class GetDistributions:
    def __init__(self, SiteGmf, EnsembleSize, Lon_Event, Lat_Event, NumGMPEsRealizations, event_dir, IMT, stationfile, imt_min, imt_max, 
                 fileScenariosWeights, pois_file, pois_subset, n_pois, max_distance, pois_selection_method, deg_round,
                 num_processes, pois_subset_flag):

        self.SiteGmf = SiteGmf
        self.EnsembleSize = EnsembleSize
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.NumGMPEsRealizations = NumGMPEsRealizations
        self.event_dir = event_dir
        self.imt = IMT
        self.stationfile = stationfile
        self.imt_min = imt_min
        self.imt_max = imt_max
        self.fileScenariosWeights = fileScenariosWeights
        self.pois_file = pois_file
        self.pois_subset = bool(pois_subset)
        self.n_pois = n_pois
        self.max_distance = max_distance
        self.pois_selection_method = pois_selection_method
        self.deg_round = deg_round
        self.num_processes = num_processes
        self.pois_subset_flag = bool(pois_subset_flag)

        if self.pois_subset_flag:
            if self.pois_subset == False:  

                self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
                self.POIs_lat = np.array(self.POIs_lat)
                self.POIs_lon = np.array(self.POIs_lon)
                print("Found ", self.n_pois, "POIs")
                self.i_pois = list(range(self.n_pois))

            else:

                if self.pois_selection_method == 'random':
                    self.azimuths = None
                    self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                                self.pois_selection_method, self.n_pois, self.max_distance)
                else:
                    self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                                self.pois_selection_method, self.n_pois, self.max_distance)
    
                self.n_pois = n_pois
                self.POIs_lat = np.array(self.POIs_lat)
                self.POIs_lon = np.array(self.POIs_lon)
                self.i_pois = self.idx_POIs

        else:
            self.pois_file = "POIs.txt"
            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths, self.n_pois = share_pois(self.pois_file)
            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)
            print("Found ", self.n_pois, "POIs in the shared file ", self.pois_file)
            self.i_pois = self.idx_POIs

    def plot_distributions(self):

        """
        Plots the CDFs at the selected POIs along with data measured at the corresponding stations
        (or at the closest stations if not available)
        """

        print("********* GETTING DISTRIBUTIONS *******")

        stats, _ = GetStatistics.calc_statistics(self)
        vector_stat = stats['vector_stat'] 
        vector = stats['vector'] 
        weight = stats['weight'] 

        data_lon, data_lat = StationRecords.get_data_coord(self)
        data_imt = StationRecords.get_data(self)
        station_id, station_name = StationRecords.get_stations(self)

        for iPoi in range(self.n_pois):
            
            sample = vector[iPoi]
            weights = weight[iPoi]
            estimator = sns.distributions.ECDF('proportion', complementary=True)
            stat, vals = estimator(sample, weights=weights)
            selPOI_CDF_vec_ecdf = stat
            selPOI_CDF_vec_vals = vals

            selPOI_p50_vec = vector_stat['Median'][iPoi]
            selPOI_p10_vec = vector_stat['Percentile 10'][iPoi]
            selPOI_p90_vec = vector_stat['Percentile 90'][iPoi]

            # SELECT THE CLOSEST RECORDING TO THE SELECTED POI
  
            selPOI_lat = self.POIs_lat[iPoi]
            selPOI_lon = self.POIs_lon[iPoi]

            tmp_dist = dist_lonlat(selPOI_lon, selPOI_lat, data_lon, data_lat, 'degree')
            
            isel = np.argmin(tmp_dist)

            datum = data_imt[isel] 
            datum_lon = data_lon[isel]
            datum_lat = data_lat[isel]
            datum_dist = tmp_dist[isel]
            station_id_isel = station_id[isel]
            station_name_isel = station_name[isel]

            print('---- POINT ', iPoi + 1, '-----')
            print('POI Coord : ', selPOI_lon,'/', selPOI_lat)
            if self.imt == 'PGV':
                print(f"Datum Coord: {datum_lon} / {datum_lat}: {self.imt} -> {datum} cm/s")
            else:
                print(f"Datum Coord: {datum_lon} / {datum_lat}: {self.imt} -> {datum} g")
           
            if datum_dist < self.max_distance: 
                datum_found = 1
            else: 
                datum_found = -1
                print('*** Observation too far (> ' + str(self.max_distance) + ' km) ***')    
            if datum < 0:
                print(f"*** No measured {self.imt} ***")

            # PLOT
            fig = plt.figure(figsize=(9,5))
            fig.patch.set_facecolor('white')

            # MAP
            plt.subplot(1,2,1)

            xlim_min = np.min(np.floor(data_lon/self.deg_round) * self.deg_round)
            xlim_max = np.max(np.ceil(data_lon/self.deg_round) * self.deg_round)
            ylim_min = np.min(np.floor(data_lat/self.deg_round) * self.deg_round)
            ylim_max = np.max(np.ceil(data_lat/self.deg_round) * self.deg_round)
            m = Basemap(projection='merc',llcrnrlat=ylim_min,urcrnrlat=ylim_max,\
                    llcrnrlon=xlim_min,urcrnrlon=xlim_max,lat_ts=20,resolution='i')
            
            x, y = m(selPOI_lon, selPOI_lat)
            m.plot(x, y, 'or', label=f"POI {iPoi + 1}")

            x, y = m(datum_lon, datum_lat)
            m.plot(x, y, 'xg', label='Closest Observation')
            
            x, y = m(self.Lon_Event, self.Lat_Event)
            m.plot(x, y, '*', label='Epicenter')

            plt.title(f"Station: {station_id_isel} ({station_name_isel})", fontsize=10)
            
            plt.legend()
            
            m.drawcoastlines()

            plt.subplot(1,2, 2)
            
            plt.plot(selPOI_CDF_vec_vals, 1-selPOI_CDF_vec_ecdf,'k-',label='CDF')
            plt.grid()
            plt.plot([datum,datum], [0,1], 'r-', label='Observation')
            if self.imt == 'PGV':
                plt.xlabel(f"{self.imt} (cm/s)")
            else:
                plt.xlabel(f"{self.imt} (g)")
            plt.ylabel("CDF")
            
            plt.semilogx([selPOI_p50_vec, selPOI_p50_vec], [0,1], 'k--', label='Median')    
            plt.semilogx([selPOI_p10_vec, selPOI_p10_vec], [0,1], 'k:', label='10th & 90th percentiles')
            plt.semilogx([selPOI_p90_vec, selPOI_p90_vec], [0,1], 'k:')
            if self.imt == 'PGV':
                plt.title(f"{self.imt} = {datum:.4f} cm/s, distance: {datum_dist:.0f} km", fontsize=10)
            else:
                plt.title(f"{self.imt} = {datum:.4f} g, distance: {datum_dist:.0f} km", fontsize=10)

            plt.legend()
            
            path = os.path.join(os.getcwd(), "OUTPUT/DISTRIBUTIONS")
            if not os.path.exists(path):
                os.makedirs(path)
            if datum > 0 and datum_found > 0:
                plt.show()
                fig.savefig(path + '/' + 'Distr_POI-' + "{:03d}".format(iPoi + 1) + '.pdf', dpi = 200)
            else:
                print('!!! Figure not saved')
                
            plt.close(fig)

        print(f"***** Figures saved in {path} *****")   

        pois_map(self.POIs_lat, self.POIs_lon, self.Lat_Event, self.Lon_Event, self.deg_round, path)


class EnsemblePlot:
    def __init__(self, SiteGmf, IMT, Lon_Event, Lat_Event, EnsembleSize, NumGMPEsRealizations, fileScenariosWeights, pois_file, 
                 pois_subset, n_pois, max_distance, deg_round, pois_selection_method, num_processes, pois_subset_flag):

        self.SiteGmf = SiteGmf
        self.EnsembleSize = EnsembleSize
        self.imt = IMT
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.NumGMPEsRealizations = NumGMPEsRealizations
        self.fileScenariosWeights = fileScenariosWeights
        self.pois_file = pois_file
        self.pois_subset = bool(pois_subset)
        self.n_pois = n_pois
        self.max_distance = max_distance     
        self.deg_round = deg_round
        self.pois_selection_method = pois_selection_method
        self.num_processes = num_processes
        self.pois_subset_flag = bool(pois_subset_flag)

        if self.pois_subset_flag:
            if self.pois_subset == False:  

                self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file)
                self.POIs_lat = np.array(self.POIs_lat)
                self.POIs_lon = np.array(self.POIs_lon)
                print("Found ", self.n_pois, "POIs")
                self.i_pois = list(range(self.n_pois))

            else:

                if self.pois_selection_method == 'random':
                    self.azimuths = None
                    self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                                self.pois_selection_method, self.n_pois, self.max_distance)
                else:
                    self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths = get_pois_subset(self.pois_file, Lon_Event, Lat_Event, 
                                                                self.pois_selection_method, self.n_pois, self.max_distance)
                    
                self.n_pois = n_pois
                self.POIs_lat = np.array(self.POIs_lat)
                self.POIs_lon = np.array(self.POIs_lon)
                self.i_pois = self.idx_POIs

        else:
            self.pois_file = "POIs.txt"
            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.azimuths, self.n_pois = share_pois(self.pois_file)
            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)
            print("Found ", self.n_pois, "POIs in the shared file ", self.pois_file)
            self.i_pois = self.idx_POIs

    def plot(self):

        """
        Returns the ensemble plot at the selected POIs
        """

        print("********* GENERATING ENSEMBLE PLOT *******")

        poi_indices = [idx + 1 for idx in range(self.n_pois)]

        stats, _ = GetStatistics.calc_statistics(self)
        vector_stat = stats['vector_stat'] 

        fig, ax = plt.subplots(figsize=(12, 6)) 
        for iPoi in range(self.n_pois):

            median = vector_stat['Median'][iPoi]
            selPOI_p5_vec = vector_stat['Percentile 5'][iPoi]
            selPOI_p95_vec = vector_stat['Percentile 95'][iPoi]
                        
            whiskerprops = dict(color='black')
            flierprops = dict(marker='o', markerfacecolor='red', markersize=8, linestyle='none')
            medianprops = dict(color='blue')
            
            ax.boxplot([[]], positions=[poi_indices[iPoi]], widths=0.6, showfliers=False, whiskerprops=whiskerprops, medianprops=medianprops, flierprops=flierprops)
            ax.vlines(poi_indices[iPoi], selPOI_p5_vec, selPOI_p95_vec, color='black')
            median_handle = ax.scatter(poi_indices[iPoi], median, color='orange', zorder=3)
            p5_handle = ax.scatter(poi_indices[iPoi], selPOI_p5_vec, color='gray', zorder=3)
            p95_handle = ax.scatter(poi_indices[iPoi], selPOI_p95_vec, color='gray', zorder=3)

        handles = [median_handle, p5_handle]
        labels = ['Median', '5th-95th percentiles']
        ax.legend(handles, labels, loc='upper right')
        if self.azimuths is None:
            ax.set_xticklabels(poi_indices)
            ax.set_xlabel('POI index', fontsize=16, labelpad=15)
        else:
            ax.set_xticklabels([int(az) for az in self.azimuths])
            ax.set_xlabel('Azimuth (°)', fontsize=16, labelpad=15)
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(poi_indices)
            ax2.set_xticklabels(poi_indices)
            ax2.set_xlabel('POI index', fontsize=16, labelpad=15)
        if self.imt == 'PGV':
            ax.set_ylabel(f"{self.imt} (cm/s)", fontsize=16, labelpad=15)  
        else:
            ax.set_ylabel(f"{self.imt} (g)", fontsize=16, labelpad=15)
   
        path = os.path.join(os.getcwd(), "OUTPUT/")
        fig.savefig(path + f"/Ensemble_Spread_Plot_{self.imt}.pdf", bbox_inches='tight')
        plt.close(fig)
        print(f"***** Figure saved in {path} *****")
 
        pois_map(self.POIs_lat, self.POIs_lon, self.Lat_Event, self.Lon_Event, self.deg_round, path)




