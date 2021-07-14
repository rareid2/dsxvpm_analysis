import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from run_rays import single_run_rays, parallel_run_rays
from bfield import getBdir
from ray_plots import plot_density_alongpath, plotrefractivesurface, plotray2D
from raytracer_utils import read_rayfile, read_input_jobs, read_damp_matlab, read_bigrayfile
from constants_settings import *
from convert_coords import convert2
import matplotlib.pyplot as plt
import random 
from ray_plots import plotray2D
import matplotlib.pyplot as plt

# get input coordinates
"""
nrays = 5013

lats, lons, psds = read_input_jobs('/home/rileyannereid/workspace/SR_output/inputs/coord_1988_fit_05.txt')
rayfile_directory = '/home/rileyannereid/workspace/SR_output' # store output here
ray_datenum = dt.datetime(2020,6,1,12,0, tzinfo=dt.timezone.utc)
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')+'_old/Kp1'

nwc_lat = 21.82

md = 6 

save_inds = []
for li, (la, lo) in enumerate(zip(lats,lons)):
    tx_loc = [R_E+500e3, la, lo] # set location

    if 8 < np.abs(tx_loc[1]) < 16:
        save_inds.append(li)
        
print('got pos')

file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []

for filename in file_titles:
    if '.ray' in filename and '05' in filename:
        this_file = filename
        print(this_file)
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

# this works!
def read_bigrayfile(rayfile):
    ''' Load output from Forest's raytracer'''
    num=0

    ray_data = []
    last_raynum = 1
    last_line = 0
    bad_ray = 0
    with open(rayfile) as a_file:

        for line in a_file:
            lines = line.split()
            ray_num = int(lines[0])

            if int(lines[1])!=1:
                if ray_num == last_raynum:
                    pass
                else:
                    bad_ray+=1

            if ray_num == last_raynum:
                pass
            else:
                ray_data.append([float(last_line[3]),float(last_line[4]),float(last_line[5])])

            last_line = lines
            last_raynum = ray_num
    print(bad_ray, ray_num)
    a_file.close()

    return ray_data

ray_data = read_bigrayfile(os.path.join(ray_out_dir, this_file))
print('read a BIG file')

coords = []
for si in save_inds:
    if si < 5012:
        coords.append(ray_data[si])

new_coords = convert2(coords, [ray_datenum for i in range(len(coords))], 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
new_lats = [nc[1] for nc in new_coords]
new_lons = [nc[2] for nc in new_coords]


#for ni,(nl,nlo) in enumerate(zip(new_lats,new_lons)):
#    if 10 < np.abs(nl) < 11:
#        print('low',ni,nl,nlo)
    #elif 37.75 < np.abs(nl) < 38:
    #    print('high',ni,nl,nlo)
    #elif 26 < np.abs(nl) < 26.5:
    #    print('mid',ni,nl,nlo)


#plt.scatter(new_lons, new_lats)
#plt.show()
#plt.close()

# high
#ray = raylist[417]
#print(len(ray['time']))
# goes from 0 to 1198

# mid
#ray = raylist[57]
#print(len(ray['time']))
# goes from 0 to 1264

# low
ray = raylist[656]
#print(len(ray['time']))
# goes from 0 to 899


for t_save in range(0,900,100):
    print(t_save)
    plotray2D(ray_datenum, [ray], ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md,t_save,show_plot=False)
    plotrefractivesurface(ray_datenum, ray,ray_out_dir,t_save)
    plot_density_alongpath(ray_datenum, ray,ray_out_dir,t_save)
"""

# get input coordinates

nworkers = 16
nrays = 20000

lats = []
lons = []
# get rand lats and lons
for n in range(nrays):
    rla  = random.random()
    # convert to latitude
    rla = (rla * 180) - 90

    rlo  = random.random()
    # convert to longitude
    rlo = (rlo * 360) - 180

    if 10 < np.abs(rla) < 50:
        lats.append(rla)
        lons.append(rlo)

nrays = len(lats)

import cartopy.crs as ccrs

ax = plt.axes(projection=ccrs.Robinson())
ax.coastlines()
ax.gridlines(draw_labels=True)
plt.scatter(lons,lats,s=3)
plt.show()
plt.close()

rayfile_directory = '/home/rileyannereid/workspace/SR_output' # store output here
ray_datenum = dt.datetime(2020,6,1,12,0, tzinfo=dt.timezone.utc)

md = 1
positions = []
dirs = []
count = 0
for la, lo in zip(lats,lons):
    tx_loc = [R_E+500e3, la, lo] # set location
    count+=1
    # convert to SM for field line tracer
    tx_loc = convert2([tx_loc], [ray_datenum], 'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
    positions.append(tx_loc[0])

positions_list = [positions[int(i * (nrays/nworkers)):int((i+1)*nrays/nworkers)] for i in range(nworkers)]

freq = 19.88e3  # Hz

# same freq and starting position for all
freqs_list = [[freq for p in range(len(d))] for d in positions_list]
directions_list = [[np.zeros(3) for p in range(len(d))] for d in positions_list]

tvec = [ray_datenum for n in range(nworkers)]
directory_list = [rayfile_directory for i in range(len(tvec))]
mds = [md for i in range(len(tvec))]

parallel_run_rays(tvec, positions_list, directions_list, freqs_list, directory_list, mds)


# create empty lists to fill with ray files and damp files
raylist = []
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
file_titles = os.listdir(ray_out_dir)
mode_name = 'mode'+str(md)
for filename in file_titles:
    if '.ray' in filename and mode_name in filename:
        this_file = filename
        print(this_file)
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

flats = []
flons = []
bad_rays = 0
ray_count = 0
for r in raylist:
    ray_count+=1
    if r['stopcond'] != 1:
        bad_rays +=1
    else:    
        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])

        final_pos = new_coords[-1]
        flats.append(final_pos[1])
        flons.append(final_pos[2])
        #if final_pos[1] > 0: 
        #    print(new_coords[0])
                
print(ray_count,bad_rays)

#ax = plt.axes(projection=ccrs.PlateCarree())
#ax.coastlines()
#ax.gridlines(draw_labels=True)
#plt.scatter(flons,flats,s=3)
#plt.show()
#plt.close()


binnum = 100
binlon = np.linspace(0,360,num=binnum)
binlat = np.linspace(-90,90,num=binnum)

# wavenormals
bin_dens = np.zeros((binnum,binnum))

for eli,(elo,ela) in enumerate(zip(flons,flats)):
    for ii, (bl1, bl2) in enumerate(zip(binlon, binlon[1:])):
        for kk, (al1, al2) in enumerate(zip(binlat, binlat[1:])):
            if bl1 < elo < bl2 and al1 < ela < al2:
                bin_dens[ii,kk] += 1

ax.set_global()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.gridlines(draw_labels=True)
ax.coastlines()
#plt.contourf(binlon,binlat, np.transpose(bin_dens),levels=3, transform=ccrs.PlateCarree(), cmap='Blues')
plt.scatter(flons,flats,s=3)
plt.show()
plt.close()

