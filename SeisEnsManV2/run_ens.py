#!/usr/bin/env python

# Import system modules
import os
import ast
import sys
import utm
import threading
import configparser
import h5py
import json
#import ray
import numpy as np
from time     import gmtime
from time     import strftime
from datetime import datetime
import shutil

sys.path.append('modules/')
# Import functions from pyPTF modules
from preload             import load_PSBarInfo
from preload             import ptf_preload
from preload             import load_Scenarios_Reg
from load_event          import load_event_parameters
from mod_parser              import parse_ptf_stdin
from mod_parser              import update_cfg
from ellipsoids          import build_location_ellipsoid_objects
from ellipsoids          import build_ellipsoid_objects
from mix_utilities        import conversion_to_utm
from lambda_bsps_load      import load_lambda_BSPS
from lambda_bsps_sep       import separation_lambda_BSPS
from pre_selection         import pre_selection_of_scenarios
from short_term                 import short_term_probability_distribution
from probability_scenarios      import compute_probability_scenarios
from ensemble_sampling_RS       import compute_ensemble_sampling_RS

##### BEGIN #####
### Configuration files preparation ###
args                                              = parse_ptf_stdin()
cfg_file                                          = args.cfg
Config                                            = configparser.RawConfigParser()
Config.read(cfg_file)
Config                                            = update_cfg(cfg=Config, args=args)
PSBarInfo                                         = load_PSBarInfo(cfg=Config, args=args)
LongTermInfo, POIs, Mesh, Region_files            = ptf_preload(cfg=Config, args=args)


### Load of the event parameters ###
event_parameters = load_event_parameters(event       = args.event,
                                         event_ingv  = args.event_ingv,
                                         format      = args.event_format,
                                         routing_key = 'INT.QUAKE.CAT',
                                         args        = args,
                                         json_rabbit = None,
                                         cfg         = Config)


RS_samp_scen     = args.nb_scen
RS_samp_run      = 1
Scenarios_PS = dict()
Scenarios_BS = dict()

print('############## Initial ensemble #################')

print('Build ellipsoids objects')
ellipses = build_ellipsoid_objects(event = event_parameters,
                                   cfg   = Config,
                                   args  = args)
print('Conversion to utm')
LongTermInfo, POIs, PSBarInfo = conversion_to_utm(longTerm  = LongTermInfo,
                                                  Poi       = POIs,
                                                  event     = event_parameters,
                                                  PSBarInfo = PSBarInfo)

##########################################################
# Set separation of lambda BS-PS
print('Separation of lambda BS-PS')
lambda_bsps = load_lambda_BSPS(cfg                   = Config,
                               args                  = args,
                               event_parameters      = event_parameters,
                               LongTermInfo          = LongTermInfo)
lambda_bsps = separation_lambda_BSPS(cfg              = Config,
                                     args             = args,
                                     event_parameters = event_parameters,
                                     lambda_bsps      = lambda_bsps,
                                     LongTermInfo     = LongTermInfo,
                                     mesh             = Mesh)

##########################################################
# Pre-selection of the scenarios
print('Pre-selection of the Scenarios')
pre_selection = pre_selection_of_scenarios(cfg                = Config,
                                           args               = args,
                                           event_parameters   = event_parameters,
                                           LongTermInfo       = LongTermInfo,
                                           PSBarInfo          = PSBarInfo,
                                           ellipses           = ellipses)


##########################################################
# # COMPUTE PROB DISTR
# print('Compute short term probability distribution')
# short_term_probability  = short_term_probability_distribution(cfg                = Config,
#                                                               args               = args,
#                                                               event_parameters   = event_parameters,
#                                                               LongTermInfo       = LongTermInfo,
#                                                               PSBarInfo          = PSBarInfo,
#                                                               lambda_bsps        = lambda_bsps,
#                                                               pre_selection      = pre_selection)


##COMPUTE PROBABILITIES SCENARIOS: line 840
# print('Compute Probabilities scenarios')
# probability_scenarios = compute_probability_scenarios(cfg                = Config,
#                                                       args               = args,
#                                                       event_parameters   = event_parameters,
#                                                       LongTermInfo       = LongTermInfo,
#                                                       PSBarInfo          = PSBarInfo,
#                                                       lambda_bsps        = lambda_bsps,
#                                                       pre_selection      = pre_selection,
#                                                       regions            = Region_files,
#                                                       Scenarios_PS       = Scenarios_PS)


############### Real sampling ########################

print('############## Creation of the ensemble  #################')

sampled_ensemble_RS = compute_ensemble_sampling_RS(cfg                = Config,
                                          args               = args,
                                          event_parameters   = event_parameters,
                                          LongTermInfo       = LongTermInfo,
                                          PSBarInfo          = PSBarInfo,
                                          lambda_bsps        = lambda_bsps,
                                          pre_selection      = pre_selection,
                                          regions            = Region_files,
                                          Scenarios_PS       = Scenarios_PS,
                                          RS_samp_scen       = RS_samp_scen)

print('############# Writing list of scenarios #############')

for Nid in range(RS_samp_run):
    RS_samp_scen=len(sampled_ensemble_RS[Nid]['real_par_scenarios_bs'][:,0])
    par=np.zeros((11))
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    filename = f'./output/list_nb{Nid}_of_{RS_samp_scen}_scenbs_{timestamp}.txt'
    myfile = open(filename, 'w')
    for Nscen in range(RS_samp_scen):
        par[:]=sampled_ensemble_RS[Nid]['real_par_scenarios_bs'][Nscen,:]
        myfile.write("%f %f %f %f %f %f %f %f %f %f\n"%(par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10]))
    myfile.close()

print('############# Copying list of scenarios to Prob Shakemap #############')

start_folder = os.path.abspath(os.path.join(os.getcwd(), './output/'))
destination_folder = os.path.abspath(os.path.join(os.getcwd(), '../INPUT_FILES/ENSEMBLE'))
files = os.listdir(start_folder)
files = [file for file in files if file.startswith('list_nb') and file.endswith('.txt')]

# Sort the files by timestamp
files.sort(key=lambda x: int(x.split('.')[0].split('_')[-1]), reverse=True)
last_file = files[0]

# Create a BACKUP of the destination folder and move all the pre-existing files there
scenarios_files = os.listdir(destination_folder)

timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
backup_folder = os.path.abspath(os.path.join(os.getcwd(), f'../INPUT_FILES/BACKUP/ENSEMBLE_{timestamp}'))
if not os.path.exists(backup_folder):
    os.makedirs(backup_folder)

for file in scenarios_files:
    destination_path = os.path.join(destination_folder, file)
    shutil.move(destination_path, os.path.join(backup_folder, file))

# Copy the new file to the destination folder
source_path = os.path.join(start_folder, last_file)
destination_path = os.path.join(destination_folder, last_file)
shutil.copy(source_path, destination_path)

print('############# DONE! #############')
