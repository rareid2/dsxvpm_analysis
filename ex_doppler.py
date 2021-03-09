import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from satellites import sat
from get_dist import antenna_MC
from bfield import getBdir
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
from constants_settings import *
from convert_coords import convert2

from shifts import dopp_delay

ray_time = dt.datetime(2020,8,17,21,20,tzinfo=dt.timezone.utc)

# change dirs to SR interface
cwd = os.getcwd()
os.chdir('/home/rileyannereid/workspace/SR_interface')

# define a satellite object
dsx = sat()            
dsx.catnmbr = 44344    
dsx.time = ray_time 
dsx.getTLE_ephem()      

vpm = sat()    
vpm.catnmbr = 45120 
vpm.time = ray_time  
vpm.getTLE_ephem()   

stf = 8200

dsx.time = ray_time
vpm.time = ray_time

vpm.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

ray_start = dsx.pos

ray_start_vel = dsx.vel[0]
ray_end_vel = vpm.vel[0]

nrays = 1e3
rayfile_directory = '/home/rileyannereid/workspace/SR_output'
# returns a vector of directions (thetas and phis must be same length) 

directions = getBdir(ray_start, ray_time, rayfile_directory, np.linspace(0,90,nrays), np.zeros(nrays))
positions = [ray_start[0] for n in range(nrays)]
freqs = [stf for n in range(nrays)]

single_run_rays(ray_time, positions, directions, freqs, rayfile_directory)

# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_time, '%Y-%m-%d %H:%M:%S')
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []
for filename in file_titles:
    if '.ray' in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

nvec_start = []
vvel_start_mag = np.sqrt(ray_start_vel[0]**2 + ray_start_vel[1]**2 + ray_start_vel[2]**2)
vvel_start = np.array([ray_start_vel[0], ray_start_vel[1], ray_start_vel[2]])/vvel_start_mag

nvec_stop = []
vvel_stop_mag = np.sqrt(ray_end_vel[0]**2 + ray_end_vel[1]**2 + ray_end_vel[2]**2)
vvel_stop = np.array([ray_end_vel[0], ray_end_vel[1], ray_end_vel[2]])/vvel_stop_mag

st_shift = []
en_shift = []

nmags = []
for r in raylist:
    rn = r['n']

    first_ind = rn.index[0]
    final_n = rn.index[-1]
    
    rtime = r['time']

    # check for bad rays
    if rtime[final_n] < 0.1: # likely did not propagate then (NEED TO CONFIRM THIS)
        continue # go to next ray

    # initial shift
    nmag = np.sqrt(rn.x[first_ind]**2 + rn.y[first_ind]**2 + rn.z[first_ind]**2)
    nmags.append(nmag)
    
    vmag = np.sqrt(ray_start_vel[0]**2 + ray_start_vel[1]**2 + ray_start_vel[2]**2)
    n_d0t_v = (rn.x[first_ind]*ray_start_vel[0] + rn.y[first_ind]*ray_start_vel[1] + rn.z[first_ind]*ray_start_vel[2])
    fshift = stf * (1 - n_d0t_v/C)
    st_shift.append(stf - fshift)

    nvec_start.append(np.array([rn.x[first_ind],rn.y[first_ind],rn.z[first_ind]])/nmag)

    # final shift
    nmag = np.sqrt(rn.x[final_n]**2 + rn.y[final_n]**2 + rn.z[final_n]**2)
    nmags.append(nmag)
    
    vmag = np.sqrt(ray_end_vel[0]**2 + ray_end_vel[1]**2 + ray_end_vel[2]**2)
    n_d0t_v = (rn.x[final_n]*ray_end_vel[0] + rn.y[final_n]*ray_end_vel[1] + rn.z[final_n]*ray_end_vel[2])
    fshift = stf * (1 - n_d0t_v/C)
    en_shift.append(stf - fshift)

    nvec_stop.append(np.array([rn.x[final_n],rn.y[final_n],rn.z[final_n]])/nmag)

    doppler_shift = fshift/1e3
    time_shifted = dt.timedelta(microseconds=rtime[final_n]) + ray_time

# get B at VPM
b0_final = getBdir(vpm.pos, ray_time, rayfile_directory, [0], [0])
b0_final = b0_final[0]
# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

ax.quiver(0,0,0,vvel_start[0],vvel_start[1],vvel_start[2],length=0.1, normalize=True,color='red',label='DSX vel')
ax.quiver(0,0,0,nvec_start[0][0],nvec_start[0][1],nvec_start[0][2],length=0.1, normalize=True,color='blue',label='B0')
ax.quiver(0,0,0,nvec_start[1][0],nvec_start[1][1],nvec_start[1][2],length=0.1, normalize=True,color='green',label='n')
ax.legend()
ax.title.set_text('Doppler Shift at DSX = '+str(round(st_shift[1],2))+ ' Hz' + '\n' + 'nmag = ' + str(round(nmags[2],2)) + '  vmag = ' + str(round(vvel_start_mag)))

ax2.quiver(0,0,0,vvel_stop[0],vvel_stop[1],vvel_stop[2],length=0.1, normalize=True,color='red',label='VPM vel')
ax2.quiver(0,0,0,b0_final[0],b0_final[1],b0_final[2],length=0.1, normalize=True,color='blue',label='B0')
ax2.quiver(0,0,0,nvec_stop[1][0],nvec_stop[1][1],nvec_stop[1][2],length=0.1, normalize=True,color='green',label='n')
ax2.legend()
ax2.title.set_text('Doppler Shift at VPM = '+str(round(en_shift[1],2))+ ' Hz' + '\n' + 'nmag = ' + str(round(nmags[3],2)) + '  vmag = ' + str(round(vvel_stop_mag)))

plt.show()
plt.close()
"""