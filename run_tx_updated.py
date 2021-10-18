# first read from the file
import json
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from ray_plots import plotray2D
from raytracer_utils import read_bigrayfile, read_rayfile
from convert_coords import convert2
import numpy as np 
import datetime as dt 
import matplotlib.pyplot as plt
R_E = 6371.2e3 

# function to find closest nodes (check ray order using initial position)
def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

# read input files from Maria
def read_input_jobs(fname):
    f = open(fname)
    lats = []
    lons = []
    psds = []
    for line in f:
        lines = line.split()
        lons.append(float(lines[0]))
        lats.append(float(lines[1]))
        psds.append(float(lines[2]))
    return lats, lons, psds

# create json config files for tx runs (need to still add in settings afterwards manually)
def create_json_file(fname, freq, alpha):
    lats, lons, psds = read_input_jobs(fname)
    # empty dictionary
    rays = []
    pos_order = []

    for la, lo in zip(lats,lons):
        ray = {"freq": freq,
            "alt": 500e3,
            "lat": la,
            "lon": lo,
            "alpha": alpha,
            "phi": 0
        }
        rays.append(ray)
        pos_order.append((la,lo))
    rays = {"rays":rays}
    return rays, pos_order

# read rayoutput
def get_ray_out(ray_out_dir):

    file_titles = os.listdir(ray_out_dir)
    mode_name = 'mode7'

    #worker_, 
    file_titles.sort(key = lambda x: int(x.split('worker_')[1][:-4]))
    raylist = []
    fcount = 0
    for filename in file_titles:
        if '.ray' in filename and mode_name in filename:
            if '_0.ray' in filename:
                pass
            else:
                print(filename)
                raylist += read_rayfile(os.path.join(ray_out_dir, filename))
                fcount+=1
            if fcount == 16:
                break
    return raylist

# get initial and final positions from ray output
def get_initial_pos(raylist):
    init_pos = []
    for r in raylist:
        ip = r['pos'].iloc[0]
        ip = [[ip.x,ip.y,ip.z]]
        pp = convert2(ip, [dt.datetime(2020,6,1,12,0,0)], 'SM','car',['m','m','m'],'GEO','sph',['m','deg','deg'])
        init_pos.append([pp[0][1],pp[0][2]])
    return init_pos

def get_final_pos(raylist):
    final_pos = []
    for r in raylist:
        ip = r['pos'].iloc[-1]
        ip = [[ip.x,ip.y,ip.z]]
        pp = convert2(ip, [dt.datetime(2020,6,1,12,0,0)], 'SM','car',['m','m','m'],'GEO','sph',['m','deg','deg'])
        final_pos.append(pp[0])
    
    return final_pos

# save output to a txt file
def write_to_txt(final_pos,fname):

    ray_count = 0
    a_file = open(fname, "w")
    for pp in final_pos:
        ray_count+=1
        final_array = [ray_count, pp[0], pp[1], pp[2]]
        line = [str(i) for i in final_array]  # convert to strings
        a_file.write(' '.join(line) + "\n")
    a_file.close()
    print('saved output to txt file for', ray_count, ' rays')

# --------------------------------------------------------------------

def output_damp(freq):
    
    file_titles = os.listdir('/home/rileyannereid/workspace/damp_'+freq)
    file_titles.sort(key = lambda x: int(x.split('worker_')[1][:-9]))

    # Using readlines()
    all_vals = []
    for filename in file_titles:
        file1 = open(os.path.join('/home/rileyannereid/workspace/damp_'+freq,filename), 'r')
        print(filename)
        new_Lines = file1.readlines()
        vals = []
        for li,line in enumerate(new_Lines):
            sline = line.split(' ')
            for sl in sline:
                if sl == '':
                    pass
                else: 
                    vals.append(sl)
        dvals = vals[1::2]
        for dv in dvals:
            dv = float(dv)
            all_vals.append(dv)

    a_file = open('/home/rileyannereid/workspace/damping_output_'+freq+'.txt', "w")
    ray_count = 0
    for dv in all_vals:
        ray_count+=1
        final_array = [ray_count, dv]
        line = [str(i) for i in final_array]  # convert to strings
        a_file.write(' '.join(line) + "\n")
    a_file.close()
# --------------------------------------------------------------------

output_damp('1988')
"""
# set input coordinates
fname = '/home/rileyannereid/workspace/dsxvpm_analysis/tx_files/coord_1988_fit_01.txt'
# set freq, and direction
freq = 19.88e3
alpha = 0
# create config file
rays, pos_order = create_json_file(fname,freq,alpha)
#with open('config.json', 'w') as json_file:
#    print('dump')
#    json.dump(rays, json_file)

# get output (after compyign from blanca)

raylist = get_ray_out('/home/rileyannereid/workspace/1988_out_newtest/2020-06-01_12_00_00')

# plot 10 rays to check them out - some of these damp to zero? 
#ray_datenum = dt.datetime(2020,6,1,12,0,0,tzinfo=dt.timezone.utc)
#for i in range(3804,3805):
#    plotray2D(ray_datenum, [raylist[i]], '/home/rileyannereid/workspace/', 'GEO', ['Re','Re','Re'], 7, show_plot=True)

# get initial and final positions
init_pos = get_initial_pos(raylist)

# confirm everything is in the correct order
for ii, (ip, po) in enumerate(zip(init_pos, pos_order)):
    latcheck = np.abs(ip[0] - po[0])
    loncheck = np.abs(ip[1] - po[1])

    if latcheck > 0.01 or loncheck > 0.01:
        if loncheck > 361.0 or loncheck < 359:
            print(latcheck,loncheck)
            print('issue at ray #',ii)

# get final position
final_pos = get_final_pos(raylist)
# save to txt file
fname = '/home/rileyannereid/workspace/output_DE_1988_new.txt'
write_to_txt(final_pos,fname)


# Using readlines()
"""
"""
file1 = open('/home/rileyannereid/workspace/output_DE_248_new.txt', 'r')
new_Lines = file1.readlines()

file2 = open('/home/rileyannereid/workspace/output_DE_248_old.txt', 'r')
old_Lines = file2.readlines()
print('files opened')
for nl, ol in zip(new_Lines, old_Lines):
    new_coords = nl.split(' ')
    old_coords = ol.split(' ')
    diff = np.abs(float(new_coords[1]) - float(old_coords[1]))
    if diff > 0.01:
        print(diff)
"""
#confirmed - 248 is in correct order as is
# so is 214
