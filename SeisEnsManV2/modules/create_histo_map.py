import sys
import h5py
import numpy as np
import pathlib
import json
import pprint
import csv
import numpy as np
import math as math
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import shutil
import datetime

def create_histo_map(**kwargs):

   Config          = kwargs.get('cfg', None)
   par             = kwargs.get('scenarios_parameters', None) 
   ee              = kwargs.get('event_parameters', None)
   args            = kwargs.get('args', None)

   ### Parameters of the event ###
   elat=ee['lat']  
   elon=ee['lon']
   edep=ee['depth']
   emag=ee['mag']
   eang=args.angles
   ename=ee['eventid']

   #Reading list
   lat=np.array(par[:,2])
   lon=np.array(par[:,3])
   minlon=np.min(lon)-1
   maxlon=np.max(lon)+1
   minlat=np.min(lat)-1
   maxlat=np.max(lat)+1
    
   ### Definition of some parameters
   parname=['magnitude', 'lat', 'lon', 'depth of the top', 'strike', 'dip', 'rake', 'area', 'length', 'average slip']
   bindata={}
   bindata['mw']=np.arange(0,10.0,0.1)
   #bindata['lon']=np.arange(np.amin(par[:,1]),np.amax(par[:,1]),0.05)
   #bindata['lat']=np.arange(np.amin(par[:,2]),np.amax(par[:,2]),0.05)
   bindata['dep']=np.arange(1,80,1)
   bindata['strike']=np.arange(0,390,15)
   bindata['dip']=np.arange(0,100,5)
   bindata['rake']=np.arange(-190,190,15)
   #bindata['area']=np.arange(np.amin(par[:,7]),np.amax(par[:,7]),20)
   #bindata['length']=np.arange(np.amin(par[:,8]),np.amax(par[:,8]),2)
   #bindata['slip']=np.arange(np.amin(par[:,9]),np.amax(par[:,9]),0.1)
   indx=['mw','dep','strike','dip','rake']
   yticks = mtick.PercentFormatter(1)
   parname=['magnitude', 'depth of the top', 'strike', 'dip', 'rake']
   lab_name=['Magnitude Mw','Depth (km)','Strike (°)','Dip (°)','Rake (°)']
   yname = ['Proportion of Mw values','Prop. of depths','Proportion of strike values','Proportion of dip values','Proportion of rake values']
    
   ### Plot for the sampled ensembles ###
   scen=len(par)
   prob=np.ones(scen)/scen
   wei=prob[:]

   fig = plt.figure(figsize=(20, 10))
    
   # Create a 2x2 grid for histograms
   ax0 = plt.subplot2grid((2, 4), (0, 0))
   ax1 = plt.subplot2grid((2, 4), (0, 1))
   ax2 = plt.subplot2grid((2, 4), (1, 0))
   ax3 = plt.subplot2grid((2, 4), (1, 1))
    
   data=par[:,1]
   ax0.hist(data,bins=bindata[indx[0]][:],weights=wei)
   ax0.yaxis.set_major_formatter(yticks)
   ax0.axvline(emag, color='red', linestyle='--', linewidth=2)
   ax0.set_ylim([0,1])
   ax0.set_xlabel(str(lab_name[0]))
   ax0.set_ylabel(str(yname[0]))
   data=par[:,5]
   ax1.hist(data,bins=bindata[indx[2]][:],weights=wei)
   ax1.yaxis.set_major_formatter(yticks)
   if eang is not None:
      ax1.axvline(eang[0], color='red', linestyle='--', linewidth=2)
   ax1.set_ylim([0,1])
   ax1.set_xlabel(str(lab_name[2]))
   ax1.set_ylabel(str(yname[2]))
   data=par[:,6]
   ax2.hist(data,bins=bindata[indx[3]][:],weights=wei)
   ax2.yaxis.set_major_formatter(yticks)
   if eang is not None:
      ax2.axvline(eang[1], color='red', linestyle='--', linewidth=2)
   ax2.set_ylim([0,1])
   ax2.set_xlabel(str(lab_name[3]))
   ax2.set_ylabel(str(yname[3]))
   data=par[:,7]
   ax3.hist(data,bins=bindata[indx[4]][:],weights=wei)
   ax3.yaxis.set_major_formatter(yticks)
   if eang is not None:
      ax3.axvline(eang[2], color='red', linestyle='--', linewidth=2)
   ax3.set_ylim([0,1])
   ax3.set_xlabel(str(lab_name[4]))
   ax3.set_ylabel(str(yname[4]))
   
   # Create a subplot for the map that spans the rows of the grid
   ax_map = plt.subplot2grid((2, 4), (0, 2), rowspan=2, colspan=2) 
   m = Basemap(llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='h',ax=ax_map)
   m.drawparallels(np.arange(-90.,91.,1.),labels=[1,0,0,0], fontsize=10)
   m.drawmeridians(np.arange(-180.,181.,1.),labels=[0,0,0,1], fontsize=10)
   m.drawcoastlines(color='k', linewidth=2)
   m.scatter(lon,lat,c='C0',s=50)
   m.scatter(elon,elat,c='C3',marker='*',s=500)

   #left, bottom, width, height = ax_map.get_position().bounds
   # Create an inset axis for the histogram overlapping the map
   #inset_ax = fig.add_axes([left + 0.1 * width, bottom + 0.7 * height, 0.3 * width, 0.2 * height])

   inset_ax = inset_axes(ax_map, width="40%", height="40%", loc='center', 
                        bbox_to_anchor=(0.58, 0.58, 0.58, 0.58), 
                        bbox_transform=ax_map.transAxes)
   
   # Plot the histogram in the inset axis
   data=par[:,4] #mag
   inset_ax.hist(data,bins=bindata[indx[1]][:],weights=wei,alpha=0.9)
   inset_ax.yaxis.set_major_formatter(yticks)
   inset_ax.axvline(edep, color='red', linestyle='--', linewidth=2)
   inset_ax.set_ylim([0,0.5])
   inset_ax.set_xlabel(str(lab_name[1]))
   inset_ax.set_ylabel(str(yname[1]))
   inset_ax.patch.set_alpha(0.5)  # Make the background transparent

   filename = './output/parameters_histo_map_'+str(ename)+'.pdf'
   plt.tight_layout()
   plt.savefig(filename,bbox_inches='tight',format='pdf')
    

def move_output():

   print('############# Copying params histos to Prob Shakemap #############')

   start_folder = os.path.abspath(os.path.join(os.getcwd(), './output/'))
   destination_folder = os.path.abspath(os.path.join(os.getcwd(), '../OUTPUT'))
   if not os.path.exists(destination_folder):
      os.makedirs(destination_folder)
   files = os.listdir(start_folder)
   filename = [file for file in files if file.startswith('parameters_histo_map_')]

   if filename:
      file = filename[0]  
   else:
      print("No parameter hysto map found")
      sys.exit()

   # Copy the new file to the destination folder
   source_path = os.path.join(start_folder, file)
   destination_path = os.path.join(destination_folder, file)
   shutil.copy(source_path, destination_path)
    
