import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from run_rays import single_run_rays, parallel_run_rays
from bfield import getBdir
from raytracer_utils import read_rayfile, read_input_jobs, read_damp_matlab, read_bigrayfile, read_bigrayfile_in
from constants_settings import *
from convert_coords import convert2
import matplotlib.pyplot as plt

from ray_plots import plotray2D

def read_omni_data(fname):
    f = open(fname)
    Kps = []
    for line in f:
        lines = line.split()
        Kp = float(lines[38])
        Kps.append(Kp)
    return Kps


# get input coordinates
nrays = 125610
lats, lons, psds = read_input_jobs('/media/rileyannereid/DEMETER/SR_output/NWC_inputs/coord_1988_fit_01.txt')
#lats = lats[100:600]
#lons = lons[100:600]

nworkers = 16
#lats = np.linspace(-41,-1,4001)
#lons = 113*np.ones(4001)

#frequencies = [ 0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. , 5.5,  6. ,  6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10. , 10.5, 11. , 11.5, 12. ]
rayfile_directory = '/home/rileyannereid/workspace' # store output here
ray_datenum = dt.datetime(2020,6,1,12,0, tzinfo=dt.timezone.utc)
# theta = 0 goes north, theta=180 goes south

#mds = [6]
#for md in mds:
#    for freq in frequencies:
#        latitudes = [-70., -63., -56., -49., -42., -35., -28., -21., -14.,  -7.,   0., 7.,  14.,  21.,  28.,  35.,  42.,  49.,  56.,  63.,  70.]
#        longitudes = [-180., -160., -140., -120., -100.,  -80.,  -60.,  -40.,  -20.,  0.,   20.,   40.,   60.,   80.,  100.,  120.,  140.,  160., 180.]
#        nrays = len(latitudes)
#        positions = []
#        dirs  = []
#        count = 0

md = 7
positions = []
pos_order = []
count = 0
for la, lo in zip(lats,lons):
    tx_loc = [R_E+500e3, la, lo] # set location
    # convert to SM car for field line tracer
    tx_loc = convert2([tx_loc], [ray_datenum], 'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
    positions.append(tx_loc[0])
    # save lat long pair
    pos_order.append((la,lo))

    #directions, ra, thetas, phis = getBdir(tx_loc, ray_datenum, rayfile_directory, thetas, np.zeros(len(thetas)), md)
    #dirs.append(directions[0])
    #print(directions[0])

    count+=1
    #print(count)

# going to have to divide by 16*3000 -- about 10.5
# nrays = 48000
# positions = positions[0:nrays]

positions_list = [positions[int(i * (nrays/nworkers)):int((i+1)*nrays/nworkers)] for i in range(nworkers)]
freq = 19.88e3  # Hz
freqs = [freq for n in range(len(positions))]
dirs = [np.zeros(3) for i in range(len(positions))]
tvec = [ray_datenum for n in range(nworkers)]

directions_list = [dirs for i in range(len(tvec))]
freqs_list = [freqs for i in range(len(tvec))]
directory_list = [rayfile_directory for i in range(len(tvec))]
mds = [md for i in range(len(tvec))]

#for p in positions_list:
#    print(len(p))
#print([(int(i * (nrays/nworkers)),int((i+1)*nrays/nworkers)) for i in range(nworkers)])

parallel_run_rays(tvec, positions_list, directions_list, freqs_list, directory_list, mds)

#single_run_rays(ray_datenum, positions, dirs, freqs, rayfile_directory, md)

# ------------------------------------- output -------------------------------------------------
# that's it! let's look at output
# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
file_titles = os.listdir(ray_out_dir)

mode_name = 'mode'+str(md)
fray_data = []
initial_raydata = []
raylist = []
x = sorted(file_titles)
for filename in x:
    if '.ray' in filename and mode_name in filename:
        #raylist += read_rayfile(os.path.join(ray_out_dir, filename))
        print(filename)
        #print(len(raylist))
        some_ray_data, init_ray_data = read_bigrayfile(os.path.join(ray_out_dir, filename))
        fray_data.append(some_ray_data)

        #init_ray_data = read_bigrayfile_in(os.path.join(ray_out_dir, filename))
        initial_raydata.append(init_ray_data)

# flatten
#ray_data = [item for sublist in fray_data for item in sublist]
#plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md)

def write_to_txt(raylist):
    
    ray_count = 0
    bad_ray = 0
    filename = ray_out_dir+'/output_mode'+str(md)+'.txt'
    a_file = open(filename, "w")
    xs = []
    ys = []
    for r in raylist:
        ray_count+=1
        #tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        #rf = tmp_coords[-1]
        new_coords = convert2([r], [ray_datenum], 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
        new_coords = new_coords[0]

        final_array = [ray_count, new_coords[0], new_coords[1], new_coords[2]]
        xs.append(new_coords[2])
        ys.append(new_coords[1])
        line = [str(i) for i in final_array]  # convert to strings
        a_file.write(' '.join(line) + "\n") 
    
    a_file.close()
    print(ray_count,bad_ray)
    colors= np.linspace(0,len(xs),num=len(xs))
    plt.scatter(xs,ys,c=colors)
    plt.show()
    plt.close()


def find_init_list(raylist):
    ray_count = 0
    final_list = []
    #FLATTEN
    flat_list = [item for sublist in raylist for item in sublist]
    #print(len(flat_list))

    for r in flat_list:
        ray_count+=1
        #tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        #rf = tmp_coords[-1]

        new_coords = convert2([r], [ray_datenum], 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
        new_coords = new_coords[0]

        final_array = (new_coords[1], new_coords[2])
        final_list.append(final_array)
    print(ray_count)
    return final_list

inital_pos = find_init_list(initial_raydata)
inital_pos = np.array(inital_pos)

# now we have the initial positions and the expected initial positions, we need to find the index that matches each
def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

flat_list_rays = [item for sublist in fray_data for item in sublist]

rays_in_order = []
#po = pos_order[424078]
# hoping that po is a tuple
# print it 
for pi,po in enumerate(pos_order):
    print(pi)
    idx = closest_node(po,inital_pos)
    rays_in_order.append(flat_list_rays[idx])
    
write_to_txt(rays_in_order)

"""
# need ray positions traced above ionosphere 1000km
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
plt.scatter(lons,lats)
plt.scatter(T_lon,T_lat)
plt.legend(['original','traced'])
plt.title('ray start points')
plt.xlabel('long [deg]')
plt.ylabel('lat [deg]')
plt.show()
plt.close()
"""

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
"""
# read single longitude output
fname = 'media/rileyannereid/DEMETER/SR_output/2020-06-01 12:00:00/output_modeGCPM_01longitude.txt'
f = open(fname)
lats = []
lons = []

for line in f:
    lines = line.split()
    lons.append(float(lines[3]))
    lats.append(float(lines[2]))
f.close()

inlats = np.linspace(-41,-1,4001)
inlons = 113*np.ones(4001)

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
plt.scatter(inlons,inlats,c='g',s=0.2)
plt.scatter(lons,lats,c='b',s=0.2)
ax.set_extent([85,135,-60,60])

ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
plt.show()
"""