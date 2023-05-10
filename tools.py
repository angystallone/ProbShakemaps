import os
import h5py
import numpy as np 
from numpy import save
import random
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np 
from scipy import constants
import json
import sys
import matplotlib.patches as mpatches
import seaborn as sns
from mpl_toolkits.basemap import Basemap
from openquake.hazardlib.geo import Point
from pyproj import Geod
geod = Geod(ellps='WGS84')


def dist_lonlat(lon1,lat1,lon2,lat2,coordtype):
    
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

    file = path + POIs_File
    print("POIs file = ", file)

    POIs_coordinates = []
    with open(file, 'r') as f:
        for line in f:
            lat = float(line.strip().split()[0])
            lon = float(line.strip().split()[1])
            POIs_coordinates.append((lat, lon))

    return POIs_coordinates  


def get_pois(POIs_File, Lon_Event, Lat_Event, POIs_Subset, pois_selection_method, n_pois, max_distance):

    path = os.path.join(os.getcwd(), "INPUT_FILES/")
    POIs_coordinates = get_pois_coordinates_from_file(path, POIs_File)
    LATs = [tup[0] for tup in POIs_coordinates]
    LONs = [tup[1] for tup in POIs_coordinates]
    print("POIs_Subset set to ", POIs_Subset)

    if POIs_Subset == False:

        POIs_lat, POIs_lon, POIs_NAMES = [], [], []

        for lat, lon in zip(LATs, LONs):

            POIs_lat.append(lat)
            POIs_lon.append(lon)
            POIs_NAMES.append(f"Site_LAT:{Point(float(lon),float(lat)).latitude}_LON:{Point(float(lon),float(lat)).longitude}")
            n_pois = len(POIs_NAMES)

        return POIs_lat, POIs_lon, POIs_NAMES, n_pois    
    
    if POIs_Subset == True:

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

            return idx_POIs, POIs_lat, POIs_lon, POIs_NAMES   

        if pois_selection_method == 'azimuth_uniform':
            DIST_BUFFER = max_distance/3
            indices = list(range(len(LONs)))

            azimuths = []
            for i in range(len(LONs)):
                azimuth, _, _ = geod.inv(LATs[i], LONs[i], Lat_Event, Lon_Event)
                azimuths.append((azimuth + 180) % 360)  # wrap azimuths to [0, 360)

            # divide the azimuth range into equal segments and select the POI closest to the center of each segment
            segment_size = 360 / n_pois
            segment_centers = np.arange(segment_size/2, 360, segment_size)

            idx_POIs, POIs_lat, POIs_lon, POIs_NAMES = [], [], [], []
            found_poi_count = 0
            for center in segment_centers:
                segment_indices = [j for j in indices if center - segment_size/2 <= azimuths[j] < center + segment_size/2]
                random.shuffle(segment_indices) 
                found_poi = False
                for j in segment_indices:
                    azimuth_diff = abs(azimuths[j] - center)
                    if azimuth_diff < segment_size/4:
                        distance = dist_lonlat(LONs[j], LATs[j], Lon_Event, Lat_Event, 'degree')
                        if max_distance - DIST_BUFFER <= distance <= max_distance + DIST_BUFFER:
                            idx_POIs.append(j)
                            POIs_lat.append(LATs[j])
                            POIs_lon.append(LONs[j])
                            POIs_NAMES.append(f"Site_LAT:{Point(float(LONs[j]), float(LATs[j])).latitude}_LON:{Point(float(LONs[j]), float(LATs[j])).longitude}")
                            found_poi = True
                            break # Exit inner loop after finding first POI
                            
                if not found_poi:
                    found_poi_count += 1
                    print(f"WARNING! No POIs have been found between {center - segment_size/2}° and {center + segment_size/2}°") 
                
            if found_poi_count != 0:
                print(found_poi_count)
                print("**** Consider decreasing max_distance value ****")
                sys.exit()       
 
            print(len(idx_POIs))    

            return idx_POIs, POIs_lat, POIs_lon, POIs_NAMES    


def get_poi_dict(POIs_NAMES):

    # Initialize a dictionary to map each POI to its index 
    poi_idx_dict = {poi: idx for idx, poi in enumerate(POIs_NAMES)}
    poi_indices = [poi_idx_dict[poi] + 1 for poi in poi_idx_dict]

    return poi_idx_dict, poi_indices


def pois_map(POIs_lat, POIs_lon, POIs_NAMES, Lat_Event, Lon_Event, deg_round, path):

    _, poi_indices = get_poi_dict(POIs_NAMES)

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
        plt.text(x[i] + 150, y[i] + 150, label, fontsize=8, color='black', weight='bold')
    m.scatter(x, y, s=12, marker='o', color='red', label="POIs")

    m.scatter(x_event, y_event, s=70, marker='*', color='blue', label="Epicenter")
    plt.legend(loc='lower left')

    plt.savefig(path + "/POIs_Map.pdf", dpi=300, bbox_inches='tight')

    plt.show()


####################################################################################
############################ PROBSHAKEMAP TOOLS ####################################
####################################################################################

class StationRecords:
    def __init__(self, event_dir, IMT, imt_min, imt_max, deg_round, stationfile):
                 
        self.event_dir = event_dir
        self.imt = IMT
        self.imt_min = imt_min
        self.imt_max = imt_max
        self.deg_round = deg_round 
        self.stationfile = stationfile

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
                tmp = data['features'][i]['properties'][key]
                if tmp != 'null':
                    if self.imt == 'PGV':
                        data_imt[i] = tmp 
                    else:
                        data_imt[i] = tmp / 100 
                else:
                    data_imt[i] = -0.01
                        
        print(f"{self.imt} range: {np.min(data_imt)} - {np.max(data_imt)}")

        return data_imt

    def plot(self):

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
        sc = m.scatter(x, y, c=np.log10(data_imt), vmin=np.log10(self.imt_min), vmax=np.log10(self.imt_max), s=5, edgecolors='black', linewidths=0.2, cmap = cm)
        cbar = plt.colorbar(sc)
        cbar.set_label(f"log10({self.imt})")

        plt.show()

        path = os.path.join(os.getcwd(), "OUTPUT")
        figname = f"/Data_stationfile_{self.imt}.pdf"
        fig.savefig(path + figname, dpi=200)   
        print("DONE!") 
        print(f"***** Figure saved in {path} *****")


class QueryHDF5:
    def __init__(self, output, Lon_Event, Lat_Event, scenario, pois_file, pois_subset, n_pois, max_distance, pois_selection_method):

        self.output = output
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.scenario = scenario
        self.pois_file = pois_file
        self.pois_subset = bool(pois_subset)
        self.n_pois = n_pois
        self.max_distance = max_distance    
        self.pois_selection_method = pois_selection_method

        if self.pois_subset == False:  

            self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            print("Found ", self.n_pois, "POIs")

        else:

            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.n_pois = n_pois

            print("Extracted ", self.n_pois, "POIs")
        

    def get_scenarios(self):   
        outputfile = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/", self.output)
        with h5py.File(outputfile, "r") as f: 
            Scenarios = list(f)
            Scenarios_n = len(Scenarios)  

        return Scenarios, Scenarios_n
    
    def get_num_gmf(self):
        outputfile = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/", self.output)
        with h5py.File(outputfile, "r") as f: 
            scen = 'Scenario_1' 
            poi = self.POIs_NAMES[0]
            GMF_n = len(list(f[scen][poi]))

        return GMF_n    

    def print_info(self):

        print('Ensemble file : ', self.output)
        print()

        _, Scenarios_n = QueryHDF5.get_scenarios(self)

        outputfile = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/", self.output)
        print(outputfile)
        with h5py.File(outputfile, "r") as f: 
            scen = 'Scenario_%d' % self.scenario
            
            GMF_n = QueryHDF5.get_num_gmf(self)

            print(f" Number of POIs: {len(self.POIs_lat)} \n Number of Scenarios: {Scenarios_n} \n Number of GMF per POI: {GMF_n}")
            
            data = {}
            for poi in self.POIs_NAMES:
                data[poi] = list(f[scen][poi])
                    
            filename = os.path.join(os.getcwd(), "OUTPUT/", 'GMF_info.txt')
            with open(filename, 'w') as f:
                for poi in self.POIs_NAMES:
                    output_str = f"GMF realizations at {poi} for {scen}: {data[poi]}\n"
                    f.write(output_str)
            print(f"GMF info for {scen} printed in the file '{filename}")

class GetStatistics:
    def __init__(self, output, Lon_Event, Lat_Event, event_dir, IMT, stationfile, imt_min, imt_max, fileScenariosWeights, pois_file, pois_subset, n_pois, max_distance, pois_selection_method, deg_round):

        self.output = output
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
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

        if self.pois_subset == False:  

            self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

            print("Found ", self.n_pois, "POIs")

        else:

            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.n_pois = n_pois

            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

            print("Extracted ", self.n_pois, "POIs")


    def calc_statistics(self):

        thresholds_stat_names = ['Mean','Mean','Median','Percentile 10','Percentile 20','Percentile 80','Percentile 90']
        thresholds_stat = np.zeros([len(thresholds_stat_names), self.n_pois])
        vector_stat = np.zeros([len(thresholds_stat_names), self.n_pois])
        thresholds_n = 6000
        thresholds_lim = np.linspace(0, 2, thresholds_n + 1)
        thresholds_cent = 0.5 * (thresholds_lim[0 : len(thresholds_lim) - 1] + thresholds_lim[1 : len(thresholds_lim)])

        Scenarios, Scenarios_n = QueryHDF5.get_scenarios(self)

        weights = np.ones(Scenarios_n) / Scenarios_n    
        if self.fileScenariosWeights != "":
            print('Weights file : ', self.fileScenariosWeights)
            with open(self.fileScenariosWeights) as f:
                lines = f.readlines()
                weights_n = len(list(lines))

                weights = np.zeros(weights_n)
                for i, line in enumerate(lines):
                    tmp = line.split()
                    weights[i] = float(tmp[0])
                f.close

        weights = weights / np.sum(weights)

        GMF_n = QueryHDF5.get_num_gmf(self)

        outputfile = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/", self.output)
        with h5py.File(outputfile, "r") as f:
            
            thresholds_distrib = np.zeros([self.n_pois, thresholds_n])   
            vector = np.zeros([self.n_pois, GMF_n * Scenarios_n])
            weight = np.zeros([self.n_pois, GMF_n * Scenarios_n]) 

            for jp, poi in enumerate(self.POIs_NAMES):

                # STEP 1: evaluate distributions mixing events in the ensemble
    
                for isc, scen in enumerate(Scenarios):

                    imt_values = np.asarray(list(f[scen][poi]))
                    
                    a, b = np.histogram(imt_values, thresholds_lim) 
                    prob = a / GMF_n
                    thresholds_distrib[jp] += prob *  weights[isc]
                    
                    ini = GMF_n * isc
                    ifi = GMF_n * (isc+1)
                    vector[jp][ini:ifi] = imt_values
                    weight[jp][ini:ifi] = weights[isc]/GMF_n 
                                        
                ## STEP 2: evaluate the statistics of the obtained distribution
                
                nsamp = 10000
                tmp = thresholds_distrib[jp]
                tmpcum = np.cumsum(tmp)
                th = np.random.rand(nsamp)
                pga_samp = np.zeros(nsamp)
                for j,samp in enumerate(th): 
                    itemindex = np.where(tmpcum > samp)
                    if np.size(itemindex) == 0:
                        itemindex_first = 0
                    else:
                        itemindex_first = itemindex[0][0]
                    pga_samp[j]= thresholds_cent[itemindex_first]

                thresholds_stat[0][jp] = np.sum(np.multiply(tmp, thresholds_cent))
                thresholds_stat[1][jp] = np.mean(pga_samp)
                thresholds_stat[2][jp] = np.median(pga_samp)
                thresholds_stat[3][jp] = np.percentile(pga_samp,10)
                thresholds_stat[4][jp] = np.percentile(pga_samp,20)
                thresholds_stat[5][jp] = np.percentile(pga_samp,80)
                thresholds_stat[6][jp] = np.percentile(pga_samp,90)
                
                vector_stat[0][jp] = np.mean(vector[jp])
                vector_stat[1][jp] = np.sum(vector[jp] * weight[jp])/np.sum(weight[jp])
                vector_stat[2][jp] = weighted_percentile(vector[jp],weight[jp],0.5)
                vector_stat[3][jp] = weighted_percentile(vector[jp],weight[jp],0.1)
                vector_stat[4][jp] = weighted_percentile(vector[jp],weight[jp],0.2)
                vector_stat[5][jp] = weighted_percentile(vector[jp],weight[jp],0.8)
                vector_stat[6][jp] = weighted_percentile(vector[jp],weight[jp],0.9) 

            stats = {} 
            stats['thresholds_stat'] = thresholds_stat
            stats['thresholds_distrib'] = thresholds_distrib
            stats['vector_stat'] = vector_stat
            stats['vector'] = vector
            stats['weight'] = weight
            stats['thresholds_cent'] = thresholds_cent

        return stats, thresholds_stat_names
    
    def save_statistics(self):

        stats, _ = GetStatistics.calc_statistics(self)
        thresholds_stat = stats['thresholds_stat'] 
        thresholds_distrib = stats['thresholds_distrib'] 
        vector_stat = stats['vector_stat'] 
        vector = stats['vector'] 
        weight = stats['weight'] 

        path = os.path.join(os.getcwd(), "OUTPUT/npyFiles")
        if not os.path.exists(path):
            os.makedirs(path)
        save(path + '/' + f"{self.imt}" + '_thresholds_stat.npy', thresholds_stat)
        save(path + '/' + f"{self.imt}" + '_thresholds_distrib.npy', thresholds_distrib)
        save(path + '/' + f"{self.imt}" + '_vector_stat.npy', vector_stat)
        save(path + '/' + f"{self.imt}" + '_vector.npy', vector)
        save(path + '/' + f"{self.imt}" + '_weight.npy', weight) 

    def plot_statistics(self):

        path = os.path.join(os.getcwd(), "OUTPUT/STATISTICS")
        if not os.path.exists(path):
            os.makedirs(path)

        stats, stat_names = GetStatistics.calc_statistics(self)
        vector_stat = stats['vector_stat'] 

        dim_point = 10

        # TARGET POINTS
        xlim_min = np.min(np.floor(self.POIs_lon/self.deg_round) * self.deg_round)
        xlim_max = np.max(np.ceil(self.POIs_lon/self.deg_round) * self.deg_round)
        ylim_min = np.min(np.floor(self.POIs_lat/self.deg_round) * self.deg_round)
        ylim_max = np.max(np.ceil(self.POIs_lat/self.deg_round) * self.deg_round)

        for j,name in enumerate(stat_names):

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
            
            tmp = vector_stat[j]
            sc = m.scatter(x, y, c=np.log10(tmp), vmin=np.log10(self.imt_min), vmax=np.log10(self.imt_max), s=dim_point, cmap = cm)
            plt.title(name)
            cbar = plt.colorbar(sc)
            cbar.set_label(f"log10({self.imt})")
            plt.show()

            figname = path + '/' + 'Stat_' + name.strip() + '.pdf'
            fig.savefig(figname, dpi = 200)
        print(f"***** Figures saved in {path} *****")


class GetDistributions:
    def __init__(self, output, Lon_Event, Lat_Event, event_dir, IMT, stationfile, imt_min, imt_max, fileScenariosWeights, pois_file, pois_subset, n_pois, max_distance, pois_selection_method, deg_round):

        self.output = output
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
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

        if self.pois_subset == False:  

            self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

            print("Found ", self.n_pois, "POIs")

        else:

            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.n_pois = n_pois

            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

            print("Extracted ", self.n_pois, "POIs")

    def plot_distributions(self):

        stats, _ = GetStatistics.calc_statistics(self)
        vector_stat = stats['vector_stat'] 
        vector = stats['vector'] 
        weight = stats['weight'] 
        thresholds_cent = stats['thresholds_cent']

        data_lon, data_lat = StationRecords.get_data_coord(self)
        data_imt = StationRecords.get_data(self)

        ipois = np.arange(0, self.n_pois, 1) 

        for iPoi in ipois:
            
            sample = vector[iPoi]
            weights = weight[iPoi]
            estimator = sns.distributions.ECDF('proportion', complementary=True)
            stat, vals = estimator(sample, weights=weights)
            selPOI_CDF_vec_ecdf = stat
            selPOI_CDF_vec_vals = vals
            
            selPOI_p50_vec = vector_stat[2][iPoi]
            selPOI_p10_vec = vector_stat[3][iPoi]
            selPOI_p90_vec = vector_stat[6][iPoi]

            # SELECT THE CLOSEST RECORDING TO THE SELECTED POI
  
            selPOI_lat = self.POIs_lat[iPoi]
            selPOI_lon = self.POIs_lon[iPoi]

            tmp_dist = dist_lonlat(selPOI_lon, selPOI_lat, data_lon, data_lat, 'degree')
            
            isel = np.argmin(tmp_dist)

            datum = data_imt[isel] 
            datum_lon = data_lon[isel]
            datum_lat = data_lat[isel]
            datum_dist = tmp_dist[isel]

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
            m.plot(x, y, 'or', label='POI')

            x, y = m(datum_lon, datum_lat)
            m.plot(x, y, 'xg', label='Closest Observation')
            
            x, y = m(self.Lon_Event, self.Lat_Event)
            m.plot(x, y, '*', label='Epicenter')
            
            plt.legend()
            
            m.drawcoastlines()

            plt.subplot(1,2, 2)
            
            plt.plot(selPOI_CDF_vec_vals, 1-selPOI_CDF_vec_ecdf,'k-',label='CDF')
            plt.grid()
            plt.plot([datum,datum], [0,1], 'r-', label='Observation')
            plt.xlabel(f"{self.imt}")
            plt.ylabel("CDF")
            
            plt.semilogx([selPOI_p50_vec, selPOI_p50_vec], [0,1], 'k--', label='Median')    
            plt.semilogx([selPOI_p10_vec, selPOI_p10_vec], [0,1], 'k:', label='10th & 90th percentiles')
            plt.semilogx([selPOI_p90_vec, selPOI_p90_vec], [0,1], 'k:')
            if self.imt == 'PGV':
                plt.title(f"{self.imt} = {datum:.4f} cm/s, distance: {datum_dist:.0f} km")
            else:
                plt.title(f"{self.imt} = {datum:.4f} g, distance: {datum_dist:.0f} km")

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

        # Create POIs map and save it
        pois_map(self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.Lat_Event, self.Lon_Event, self.deg_round, path)


class EnsemblePlot:
    def __init__(self, output, IMT, Lon_Event, Lat_Event, pois_file, pois_subset, n_pois, max_distance, deg_round, pois_selection_method):

        self.output = output
        self.imt = IMT
        self.Lon_Event = Lon_Event
        self.Lat_Event = Lat_Event
        self.pois_file = pois_file
        self.pois_subset = bool(pois_subset)
        self.n_pois = n_pois
        self.max_distance = max_distance     
        self.deg_round = deg_round
        self.pois_selection_method = pois_selection_method

        if self.pois_subset == False:  

            self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.n_pois = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

            print("Found ", self.n_pois, "POIs")

        else:

            self.idx_POIs, self.POIs_lat, self.POIs_lon, self.POIs_NAMES = get_pois(self.pois_file, self.Lon_Event, self.Lat_Event, self.pois_subset, 
                                                              self.pois_selection_method, self.n_pois, self.max_distance)

            self.n_pois = n_pois

            self.POIs_lat = np.array(self.POIs_lat)
            self.POIs_lon = np.array(self.POIs_lon)

            print("Extracted ", self.n_pois, "POIs")


    def plot(self):

        path = os.path.join(os.getcwd(), "OUTPUT/")
        # Create POIs map and save it
        pois_map(self.POIs_lat, self.POIs_lon, self.POIs_NAMES, self.Lat_Event, self.Lon_Event, self.deg_round, path)

        _, EnsembleSize = QueryHDF5.get_scenarios(self)

        poi_idx_dict, poi_indices = get_poi_dict(self.POIs_NAMES)

        outputfile = os.path.join(os.getcwd(), "OUTPUT/HDF5_FILES/", self.output)
        with h5py.File(outputfile, "r") as f:
            # Loop over all scens
            mean_all_scens = np.zeros((self.n_pois, EnsembleSize))
            for i, scen in enumerate(f.keys()):
                # Initialize an empty array to store the mean PGA values for each POI for the current scenario
                mean_scen = np.zeros((self.n_pois,))
                #print(self.n_pois)

                # Loop over POIs in the current scenario
                for poi in f[scen]:
                    if poi in poi_idx_dict:
                        # Get the distribution of IMT values for the POI if in list POIs
                        pga_dist = np.array(f[scen][poi])
                        # Calculate the mean IMT value for the current POI
                        mean_scen[poi_idx_dict[poi]] += np.mean(pga_dist)

                # Store the mean IMT values for all POIs for the current scenario
                mean_all_scens[:, i] = mean_scen        
    
            # Create the ensemble spread plot
            fig, ax = plt.subplots()
            bp = ax.boxplot(mean_all_scens.T, positions=poi_indices, sym='', whis=[5, 95], widths=0.6)
            # set x-axis labels to POI indexes
            ax.set_xticklabels(poi_indices)

            median_patch = mpatches.Patch(color='orange', label='Median')
            whisker_patch = mpatches.Patch(color='gray', label='5th-95th Percentiles')
            plt.legend(handles=[median_patch, whisker_patch], loc='upper right')

            # set axis labels and title
            ax.set_xlabel('POI index', fontsize=16)
            if self.imt == 'PGA':
                ax.set_ylabel(f"{self.imt} mean (g)", fontsize=16)
            if self.imt == 'PGV':
                ax.set_ylabel(f"{self.imt} mean (cm/s)", fontsize=16)

            # show plot
            plt.show()
            fig.savefig(path + f"/Ensemble_Spread_Plot_{self.imt}.pdf", bbox_inches='tight')
            print(f"***** Figure saved in {path} *****")
 
        


