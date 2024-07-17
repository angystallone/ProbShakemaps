import requests
import xml.etree.ElementTree as ET
import os
from datetime import datetime, timedelta
import json
import numpy as np
from obspy import read_events
from obspy.core.utcdatetime import UTCDateTime
import shutil
import sys

# LOAD INFO FROM event.xml FILE

input_dir = os.path.abspath(os.path.join(os.getcwd(), '../INPUT_FILES'))
data_dir = os.path.abspath(os.path.join(os.getcwd(), '../data'))
input_file = os.path.join(input_dir, 'input_file.txt')

with open(input_file, 'r') as file:
    lines = file.readlines()    
    ID_Event = str(lines[3].split( )[1])

event_dir = os.path.join(data_dir, ID_Event, "current")
if not os.path.exists(event_dir):
    raise NotADirectoryError(f"{event_dir} is not a valid directory.")

eventxml = os.path.join(event_dir, "event.xml")    
if not os.path.isfile(eventxml):
    raise FileNotFoundError(f"{eventxml} does not exist.")   

tree = ET.parse(eventxml)
root = tree.getroot()
mag=root.attrib.get('mag') 
time=root.attrib.get('time')[:-1]
timestamp = datetime.strptime(time, "%Y-%m-%dT%H:%M:%S")

starttime = timestamp - timedelta(minutes=30)
endtime = timestamp + timedelta(minutes=30)
min_mag = float(mag) - 0.1

# DOWNLOAD QUAKEML FILE

url = f"https://webservices.ingv.it/fdsnws/event/1/query?starttime={starttime}&endtime={endtime}&format=quakeml&minmag={min_mag}"
print(url)
eventquakeml = 'event.quakeml'

try:
    response = requests.get(url)
    if response.status_code == 200:
        with open(eventquakeml, 'wb') as file:
            file.write(response.content)
        print(f'############# Event info downloaded to {eventquakeml} #############')
    else:
        print(f'Failed to download. Status code: {response.status_code}')
        sys.exit(1) 
except requests.exceptions.RequestException as e:
    print(f'An error occurred: {e}')
    sys.exit(1) 

# LOAD NEEDED INFO FROM .quakeml file

try:
    catalog = read_events(eventquakeml)

    # Check if catalog contains any events
    if catalog:
        event = catalog[0]  
        description = event.event_descriptions[0]  
        region_name = description.text

        try:
            origin = event.origins[0]
            time = origin.time
            utc_time = UTCDateTime(time)
            time_str = utc_time.strftime("%Y-%m-%dT%H:%M")

            latitude = origin.latitude
            longitude = origin.longitude
            depth = origin.depth / 1000.0  # Convert depth to kilometers

            magnitude = event.magnitudes[0]
            mag_type = magnitude.magnitude_type
            mag = magnitude.mag

            print(f"Event found: {region_name}, {time_str}, {latitude}, {longitude}, {depth} km, {mag_type} {mag}")

        except IndexError:
            print("Error: No origin information found for the event.")
            sys.exit(1) 

    else:
        print("No events found in catalog.")

except Exception as e:
    print(f"Error occurred: {str(e)}")
    sys.exit(1)  

print("Region name:", region_name)
print("Time:", time_str)
print("Latitude:", latitude)
print("Longitude:", longitude)
print("Depth (km):", depth)    
print("Magnitude:", mag)  
print("Magnitude type:", mag_type)  

# BUILD DEFAULT COVARIANCE MATRIX 

# (units: km)
xx = 10.0
yy = 10.0
zz = 10.0
xy = 0.0
xz = 0.0
yz = 0.0

# DEFAULT MAG PERCENTILES

p50 = mag
p16 = mag - 0.3
p84 = mag + 0.3

# BACKUP of 'event_stat.json'
file_dir = os.path.join(os.getcwd(), 'input')
eventstatjson = 'event_stat.json'
eventstatjson_fullpath = os.path.join(file_dir, eventstatjson)
print("eventstatjson_fullpath: ", eventstatjson_fullpath)

timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
if os.path.exists(eventstatjson_fullpath):
    backup_folder = os.path.join(os.path.join(os.getcwd()), 'BACKUP/EVENT_STAT')
    if not os.path.exists(backup_folder):
        os.makedirs(backup_folder)

    backup_filename = f"event_stat_{timestamp}.json"
    backup_full_path = os.path.join(backup_folder, backup_filename)

    shutil.copy(eventstatjson_fullpath, backup_full_path)

# BACKUP OF 'output' FOLDER
start_folder = os.path.join(os.getcwd(), 'output')
if not os.path.exists(start_folder):
    os.makedirs(start_folder)
timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
backup_folder = os.path.join(os.path.join(os.getcwd()), f"BACKUP/LIST_SCENARIOS_{timestamp}")
if not os.path.exists(backup_folder):
    os.makedirs(backup_folder)
for item in os.listdir(start_folder):
    item_path = os.path.join(start_folder, item)
    if os.path.isfile(item_path) and item != ".DS_Store":
        shutil.move(item_path, backup_folder)

# WRITE event_stat.json FILE
with open(eventstatjson_fullpath, 'r') as file:
    data = json.load(file)

coordinates = [longitude, latitude, depth]
data['features'][0]['geometry']['coordinates'] = coordinates
data['features'][0]['properties']['time'] = time_str
data['features'][0]['properties']['mag'] = mag
data['features'][0]['properties']['magType'] = mag_type
data['features'][0]['properties']['place'] = region_name
data['features'][0]['properties']['mag_percentiles']['p16'] = p16
data['features'][0]['properties']['mag_percentiles']['p50'] = p50
data['features'][0]['properties']['mag_percentiles']['p84'] = p84
data['features'][0]['properties']['cov_matrix']['XX'] = xx
data['features'][0]['properties']['cov_matrix']['XY'] = xy
data['features'][0]['properties']['cov_matrix']['XZ'] = xz
data['features'][0]['properties']['cov_matrix']['YY'] = yy
data['features'][0]['properties']['cov_matrix']['YZ'] = yz
data['features'][0]['properties']['cov_matrix']['ZZ'] = zz

with open(eventstatjson_fullpath, 'w') as file:
    json.dump(data, file, indent=4) 

print(f'############# Stats in {eventstatjson} updated #############')
print('############# Sampling list of scenarios #############')