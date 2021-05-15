import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from satellites import sat
from get_dist import antenna_MC
from bfield import getBdir, getBline
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile, read_input_jobs, read_damp_matlab
from constants_settings import *
from convert_coords import convert2
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt 
from transmitters import vlf_tx
import xlwt
from xlwt import Workbook
  

import random
from ray_plots import plotray2D, plotrefractivesurface, plot_damp, plotgeomfactor, stix_parameters, plot_plasmasphere


# get input coordinates
nrays = 346
lats, lons, psds = read_input_jobs('/home/rileyannereid/workspace/SR_output/inputs/coord_1988.txt')

rayfile_directory = '/home/rileyannereid/workspace/SR_output' # store output here

ray_datenum = dt.datetime(2020,6,1,12,0, tzinfo=dt.timezone.utc)

# need ray positions traced above ionosphere 1000km
crs='GEO'
carsph='sph'
units=['m','deg','deg'] 
positions = []
T_converts = []
save_for_later = []
for la, lo in zip(lats,lons):
    tx_loc = [R_E+475e3, la, lo] # set location
    # convert to SM car for field line tracer
    tx_loc = convert2([tx_loc], [ray_datenum], 'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
    # trace north to 1000 km
    T = getBline(tx_loc[0], ray_datenum, 1000)
    
    # repack
    T_repackx = T.x
    T_repacky = T.y
    T_repackz = T.z
    T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
    LshellT = [ray_datenum for i in range(len(T_repack))]
    
    T_convert = convert2([T_repack[-2]], [ray_datenum], 'SM','car',['Re','Re','Re'], 'SM', 'car', ['m','m','m'])
    T_save = convert2([T_repack[-2]], [ray_datenum], 'SM','car',['Re','Re','Re'], 'GEO', 'sph', ['m','deg','deg'])
    save_for_later.append(T_save[0])
    T_converts.append(T_convert[0])

    T_lat = [tt[1] for tt in T_converts]
    T_lon = [tt[2] for tt in T_converts]

    positions.append(tx_loc[0])

"""
plt.scatter(lons,lats)
plt.scatter(T_lon,T_lat)
plt.legend(['original','traced'])
plt.title('ray start points')
plt.xlabel('long [deg]')
plt.ylabel('lat [deg]')
plt.show()
plt.close()
"""

# theta = 0 goes north, theta=180 goes south
directions = [[0,0,0] for n in range(nrays)]

# only mds 1,6,7 currently working
md = 7
freq = 19.88e3  # Hz
freqs = [freq for n in range(nrays)]

#single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md)

# ------------------------------------- output -------------------------------------------------
# that's it! let's look at output
# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []
damplist = []
for filename in file_titles:
    if '.ray' in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', ['Re','Re','Re'],md)
#plot_plasmasphere(md)


def write_output_excel(raylist):
    # Workbook is created
    wb = Workbook()
    
    # add_sheet is used to create sheet.
    sheet1 = wb.add_sheet('Sheet 1')
    
    sheet1.write(0, 0, 'ray num')
    sheet1.write(0, 1, 'ray initial lat')
    sheet1.write(0, 2, 'ray initial lon')
    sheet1.write(0, 3, 'ray initial alt')
    sheet1.write(0, 4, 'ray start psd')
    sheet1.write(0, 5, 'ray stopcond')
    sheet1.write(0, 6, 'ray time')
    sheet1.write(0, 7, 'ray lat')
    sheet1.write(0, 8, 'ray lon')
    sheet1.write(0, 9, 'ray alt')

    sheet1.write(0, 10, 'landau damping factor')
    sheet1.write(0, 11, 'ray freq')

    ray_count = 0
    ii=0
    sheets = 2
    for r in raylist:
        ray_count+=1
        print(ray_count)

        # comes in as SM car in m 
        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        #tmp_kcoords = list(zip((w/C) * r['n'].x, (w/C) * r['n'].y, (w/C) * r['n'].z))
        
        # convert to a unit vector first
        #unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]

        #new_kcoords = convert2(unitk, tvec_datetime, 'SM', 'car', ['Re','Re','Re'], crs, carsph, ['Re','Re','Re'])
        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

        # save it
        for ni,nc in enumerate(new_coords):
            ii+=1
            sheet1.write(ii, 0, ray_count)
            sheet1.write(ii, 1, lats[ray_count-1])
            sheet1.write(ii, 2, lons[ray_count-1])
            sheet1.write(ii, 3, R_E+475e3)
            sheet1.write(ii, 4, psds[ray_count-1])
            sheet1.write(ii, 5, int(r['stopcond']))
            sheet1.write(ii, 6, tvec_datetime[ni].strftime("%m/%d/%Y, %H:%M:%S.%f"))
            sheet1.write(ii, 7, nc[1])
            sheet1.write(ii, 8, nc[2])
            sheet1.write(ii, 9, nc[0])
            sheet1.write(ii, 10, 0)
            sheet1.write(ii, 11, freq)

        if ii > 63000:
            # finish the ray and start the next sheet
            sheet1 = wb.add_sheet('Sheet' + str(sheets))
            # reset ii
            ii = 1
            sheets+=1


    wb.save(ray_out_dir+'/output_mode'+str(md)+'.xlsx')
    return ii

def write_to_txt(raylist):
    ray_count = 0
    filename = ray_out_dir+'/output_mode'+str(md)+'.txt'
    a_file = open(filename, "w")
    print(len(raylist))
    for r in raylist:
        ray_count+=1
        #print(ray_count)

        # comes in as SM car in m 
        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        #tmp_kcoords = list(zip((w/C) * r['n'].x, (w/C) * r['n'].y, (w/C) * r['n'].z))
        
        # convert to a unit vector first
        #unitk = [(tmp_kcoords[s][0], tmp_kcoords[s][1], tmp_kcoords[s][2]) / np.sqrt(tmp_kcoords[s][0]**2 + tmp_kcoords[s][1]**2 + tmp_kcoords[s][2]**2) for s in range(len(tmp_kcoords))]
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]

        #new_kcoords = convert2(unitk, tvec_datetime, 'SM', 'car', ['Re','Re','Re'], crs, carsph, ['Re','Re','Re'])
        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)

        # save it
        nc = new_coords[-1]
        final_array = [ray_count, nc[0], nc[1], nc[2]]
        line = [str(i) for i in final_array]  # convert to strings
        a_file.write(' '.join(line) + "\n") 


    a_file.close()


write_to_txt(raylist)

def write_damp_to_excel(ray_out_dir):
    file_titles = os.listdir(ray_out_dir)

    dd = read_damp_matlab(os.path.join(ray_out_dir,'ray_out_mode1_damp.txt'))
    
    return dd
