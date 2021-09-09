import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from run_rays import single_run_rays, parallel_run_rays
from ray_plots import plotray2D
from bfield import getBdir
from raytracer_utils import read_rayfile, read_input_jobs, read_damp_matlab, read_bigrayfile, read_bigrayfile_in
from constants_settings import *
from convert_coords import convert2
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter
import PyGeopack as gp

# script to run ray tracing from VLF transmitters
# sort the output

# ...................................................helper funcs...................................
def write_to_txt(raylist):
    
    ray_count = 0
    bad_ray = 0
    filename = ray_out_dir+'/output_mode'+str(md)+'.txt'
    a_file = open(filename, "w")
    xs = []
    ys = []
    for r in raylist:
        ray_count+=1

        new_coords = convert2([r], [ray_datenum], 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
        new_coords = new_coords[0]
        xs.append(new_coords[2])
        ys.append(new_coords[1])

        final_array = [ray_count, new_coords[0], new_coords[1], new_coords[2]]
        line = [str(i) for i in final_array]  # convert to strings
        a_file.write(' '.join(line) + "\n")

    print('saved output for', ray_count, 'rays')
    print(len(xs))
    plt.scatter(xs,ys)
    plt.show()
    plt.close() 
    
    a_file.close()


def find_init_list(raylist):
    ray_count = 0
    final_list = []
    xs = []
    ys = []
    #FLATTEN
    flat_list = [item for sublist in raylist for item in sublist]

    for r in flat_list:
        ray_count+=1

        new_coords = convert2([r], [ray_datenum], 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
        new_coords = new_coords[0]

        if new_coords[2] > 180:
            long_new = new_coords[2] - 180
            long_new = -1*(180 - long_new)
        else:
            long_new = new_coords[2]

        final_array = (float(new_coords[1]), float(long_new))
        final_list.append(final_array)
        xs.append(long_new)
        ys.append(new_coords[1])
    plt.scatter(xs,ys)
    plt.show()
    plt.close()
    print('found initial positions for', ray_count, 'rays')
    return final_list

def read_omni_data(fname):
    f = open(fname)
    Kps = []
    for line in f:
        lines = line.split()
        Kp = float(lines[38])
        Kps.append(Kp)
    return Kps

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)
# ...............................................................................................................

# get input coordinates
nrays = 70665
lats, lons, psds = read_input_jobs('/media/rileyannereid/DEMETER/SR_output/tx_inputs/coord_214_fit_01_cut.txt')

nrays = len(lats)

# plot to double check them
plt.scatter(lons,lats)
plt.show()
plt.close()

#set up
nworkers = 16
md = 7
freq = 21.4e3  # Hz
ray_datenum = dt.datetime(2020,6,1,12,0, tzinfo=dt.timezone.utc)
rayfile_directory = '/home/rileyannereid/workspace' # store output here

# save positions here
positions = []
pos_order = []
dirs = []
count = 0

makedate = ray_datenum.strftime('%Y%m%d')
Date = int(makedate)
ut = ray_datenum.hour + ray_datenum.minute/60 + ray_datenum.second/3600

for la, lo in zip(lats,lons):
    tx_loc = [R_E+500e3, la, lo] # set location

    # convert to SM car for field line tracer
    tx_loc = convert2([tx_loc], [ray_datenum], 'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
    positions.append(tx_loc[0])

    # save lat long pair
    pos_order.append((la,lo))
    count+=1

    # get Bfield direction
    Bx,By,Bz = gp.ModelField(tx_loc[0][0]/R_E,tx_loc[0][1]/R_E,tx_loc[0][2]/R_E,Date,ut,Model='T96',CoordIn='SM',CoordOut='SM')
    Bdir = np.array([Bx, By, Bz])
    Bunit = Bdir/np.linalg.norm(Bdir)
    Bsouth = [-1*float(Bunit[0]), -1*float(Bunit[1]), -1*float(Bunit[2])]
    dirs.append(Bsouth)

    # for a quick run
    #dirs.append([0,0,0])

print('starting',count, 'rays')
            
# chunk up the starting list for each worker
positions_list = [positions[int(i * (nrays/nworkers)):int((i+1)*nrays/nworkers)] for i in range(nworkers)]
directions_list = [dirs[int(i * (nrays/nworkers)):int((i+1)*nrays/nworkers)] for i in range(nworkers)]

freqs = [freq for n in range(len(positions))]
tvec = [ray_datenum for n in range(nworkers)]

# all start field aligned
#dirs = [np.zeros(3) for i in range(len(positions))]
#directions_list = [dirs for i in range(len(tvec))]

freqs_list = [freqs for i in range(len(tvec))]
directory_list = [rayfile_directory for i in range(len(tvec))]
mds = [md for i in range(len(tvec))]

# --------------------------------------------- run ---------------------------------------------
parallel_run_rays(tvec, positions_list, directions_list, freqs_list, directory_list, mds)
#single_run_rays(tvec[0], positions_list[0], directions_list[0], freqs_list[0], directory_list[0], mds[0])

# --------------------------------------------- output ------------------------------------------
# that's it! let's look at output
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
file_titles = os.listdir(ray_out_dir)
mode_name = 'mode'+str(md)

fray_data = []
initial_raydata = []
#nraylist = []
x = sorted(file_titles)
for filename in x:
    if '.ray' in filename and mode_name in filename:
        #nraylist += read_rayfile(os.path.join(ray_out_dir, filename))
        #print(filename)

        some_ray_data, init_ray_data = read_bigrayfile(os.path.join(ray_out_dir, filename))
        fray_data.append(some_ray_data)
        initial_raydata.append(init_ray_data)
#plotray2D(ray_datenum, nraylist[:200], ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md, show_plot=True,damping_vals=None)

# --------------------------------------------- sort ------------------------------------------
# get the starting lats and lons for each ray
inital_pos = find_init_list(initial_raydata)
inital_pos = np.array(inital_pos)

# create a flat list
flat_list_rays = [item for sublist in fray_data for item in sublist]

rays_in_order = [] 
# based on input initial positions and extracted initial positions, sort the rays
for pi,po in enumerate(pos_order):
    idx = closest_node(po,inital_pos)
    rays_in_order.append(flat_list_rays[idx])

# finally, plot and save
write_to_txt(rays_in_order)
