import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from satellites import sat
from bfield import getBdir, getBline
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
from constants_settings import *
from convert_coords import convert2
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt 
import random
from random import randrange, uniform
import matplotlib as mpl
from geopy.distance import geodesic
from ray_plots import plotray2D, plotrefractivesurface, plotgeomfactor, stix_parameters
from mpl_toolkits.axes_grid1 import make_axes_locatable
from michaelcode.OutputReader import OutputReader

#GetRayGroup(ray number) [indexed starting at 1!], and then do a GetRay(8200) [indexed by frequency], you should have a Ray object with tg, r and n in it.
rg = OutputReader('/home/rileyannereid/workspace/dsxvpm_analysis/michaelcode/rrays.h5')
pp = rg.GetRayGroup(1)
pg = pp.GetRay(8200)
#print(pg.r)

md = 7
nrays = 10
nworkers = 1

freq = 8.2e3
rayfile_directory = '/media/rileyannereid/DEMETER/SR_output' # store output here

ray_datenum = dt.datetime(2020,5,19,15,47,45,tzinfo=dt.timezone.utc)
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')

# ray directions to use
d0 = [7.070570306838824e-01, 3.662883097385405e-01, -6.048993548594463e-01] 
d1 = [7.128028910778679e-01, 7.004703186830318e-01, 3.540298172647105e-02]
d2 = [8.854369009445802e-01, -2.491765283955410e-01, 3.923169026978113e-01]
d3 = [5.070268795884871e-01, 5.060436221164279e-01, 6.977417831046437e-01]
d4 = [-2.228416860601109e-02, -8.929300558652086e-01, 4.496435601252332e-01]
d5 = [3.687740382211073e-01, 9.106685786106854e-02, 9.250473156187250e-01]
d6 = [6.008230239445215e-01, -7.391276285655751e-01, 3.044701013058425e-01]
d7 = [7.805500877221502e-01, -2.615104194281180e-01, 5.677621518624446e-01]
d8 = [9.496430544885617e-01, 1.462538061920225e-01, 2.771062850892790e-01]
d9 = [8.482729000399917e-01, -2.021793691901503e-01, 4.894451856251273e-01]

# we need the positions of the satellites -- use the sat class
dsx = sat()             # define a satellite object
dsx.catnmbr = 44344     # provide NORAD ID
dsx.time = ray_datenum  # set time
dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

# propagate the orbit! setting sec=0 will give you just the position at that time
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='GEO', carsph='car', units=['km','km','km'])

# ray start position
ray_start = dsx.pos

dz = [0,0,0]
directions_list = [[d0,d1,d2,d3,d4,d5,d6,d7,d8,d9]]

directions = [dz for i in range(10)]

# same freq and starting position for all
freqs_list = [[freq for p in range(len(d))] for d in directions_list]
positions = [[ray_start[0] for p in range(len(d))] for d in directions_list]

positions = [ray_start[0] for p in range(10)]
freqs = [freq for p in range(10)]

tvec = [ray_datenum for n in range(nworkers)]
directory_list = [rayfile_directory for i in range(len(tvec))]
mds = [md for i in range(len(tvec))]
print(positions)
#parallel_run_rays(tvec, positions, directions_list, freqs_list, directory_list, mds)
single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md, runmodeldump=False)

raylist = []
mode_name = 'mode' + str(md)
# use mode name to avoid workers of the same label
file_titles = os.listdir(ray_out_dir)

x = sorted(file_titles)
for filename in x:
    if '.ray' in filename and mode_name in filename:
        if 'Main' in filename: # avoid call from bfield
            pass
        else:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))
            print(filename)

# plots relative power 2D, returns maximum refractive index from ra (uses Stix params function)

for i in range(10):
    nmax, max_wna = plotray2D(ray_datenum, [raylist[i]], ray_out_dir, 'MAG', ['Re','Re','Re'], md, show_plot=True,return_nmax=False, checklat=None)
#ray_datenum, raylist, ray_out_dir, crs, units, md, show_plot=True, return_nmax=False, checklat=None

