import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from satellites import sat
from get_dist import antenna_MC
from bfield import getBdir
from bfield import getBline

from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
from constants_settings import *
from convert_coords import convert2

import random

# function to find doppler shifts from satellite velocities


def get_lshell_sats(stime, ss):
    os.chdir('/home/rileyannereid/workspace/SR_interface')
    stime = dt.datetime(stime.year,stime.month,stime.day,stime.hour,stime.minute,stime.second,tzinfo=dt.timezone.utc)

    mysat = sat()             # define a satellite object
    mysat.catnmbr = ss        # provide NORAD ID
    mysat.time = stime        # set time
    mysat.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

    # propagate the orbit! setting sec=0 will give you just the position at that time
    mysat.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

    satpos = mysat.pos
    thatline = getBline(satpos[0], stime)

    return thatline.Lshell

def dopp_delay(nrays, rayfile_directory, tnt_times_shift, dur_shift, startf_shift, stopf_shift):

    find_tle_time = tnt_times_shift[0]
    find_tle_time = find_tle_time.replace(tzinfo=dt.timezone.utc)

    # get angle defs
    thetas, phis = antenna_MC(nrays)

    thetas = []
    for nr in range(nrays//2):
        th = random.randrange(60, 90)
        th = random.randrange(-90, -60)
        thetas.append(th)
    
    phis = np.zeros(nrays)

    # change dirs to SR interface
    cwd = os.getcwd()
    os.chdir('/home/rileyannereid/workspace/SR_interface')
    
    # define a satellite object
    dsx = sat()            
    dsx.catnmbr = 44344    
    dsx.time = find_tle_time 
    dsx.getTLE_ephem()      

    vpm = sat()    
    vpm.catnmbr = 45120 
    vpm.time = find_tle_time  
    vpm.getTLE_ephem()    

    # loop through tnt times 
    tnt_dop = []
    tnt_t = []

    #record all shifts
    alldop = []
    allsec = []
    allthetas = []
    for tim, dur, strf, stpf in zip(tnt_times_shift, dur_shift, startf_shift, stopf_shift):

        pulse_t = []
        pulse_freqs = []

        if dur == 150:
            pulse_t.append(tim)
            pulse_t.append(tim + dt.timedelta(microseconds=75e3))
            pulse_t.append(tim + dt.timedelta(microseconds=150e3))

            pulse_freqs.append(strf)
            pulse_freqs.append(strf+100)
            pulse_freqs.append(strf-100)

        elif dur == 250: # exclude large f ramps
            pulse_t.append(tim)
            pulse_t.append(tim+dt.timedelta(microseconds=dur*1e3))
            pulse_freqs.append(strf)
            pulse_freqs.append(stpf)
        
        # loop through 'pulses'
        pulse_dop = []
        pulse_tdelay = []
        for t_time, freq in zip(pulse_t, pulse_freqs):

            dsx.time = t_time
            vpm.time = t_time

            vpm.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])
            dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

            ray_start = dsx.pos

            ray_start_vel = dsx.vel[0]
            ray_end_vel = vpm.vel[0]

            # returns a vector of directions (thetas and phis must be same length) 
            directions = getBdir(ray_start, t_time, rayfile_directory, thetas, phis)
            positions = [ray_start[0] for n in range(nrays)]
            freqs = [freq for n in range(nrays)]

            single_run_rays(t_time, positions, directions, freqs, rayfile_directory)

            # Load all the rayfiles in the output directory
            ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(t_time, '%Y-%m-%d %H:%M:%S')
            file_titles = os.listdir(ray_out_dir)

            # create empty lists to fill with ray files and damp files
            raylist = []
            for filename in file_titles:
                if '.ray' in filename:
                    raylist += read_rayfile(os.path.join(ray_out_dir, filename))

            doppler_shifted = []
            time_shifted = []
            for ri, r in enumerate(raylist):
                rn = r['n']

                first_ind = rn.index[0]
                final_n = rn.index[-1]

                rtime = r['time']

                # check for bad rays
                if rtime[final_n] < 0.01: # likely did not propagate then (NEED TO CONFIRM THIS)
                    continue # go to next ray

                # initial shift
                nmag = np.sqrt(rn.x[first_ind]**2 + rn.y[first_ind]**2 + rn.z[first_ind]**2)
                vmag = np.sqrt(ray_start_vel[0]**2 + ray_start_vel[1]**2 + ray_start_vel[2]**2)
                n_d0t_v = (rn.x[first_ind]*ray_start_vel[0] + rn.y[first_ind]*ray_start_vel[1] + rn.z[first_ind]*ray_start_vel[2])
                fshift = freq * (1 - n_d0t_v/C)

                # final shift
                nmag = np.sqrt(rn.x[final_n]**2 + rn.y[final_n]**2 + rn.z[final_n]**2)
                vmag = np.sqrt(ray_end_vel[0]**2 + ray_end_vel[1]**2 + ray_end_vel[2]**2)
                n_d0t_v = (rn.x[final_n]*ray_end_vel[0] + rn.y[final_n]*ray_end_vel[1] + rn.z[final_n]*ray_end_vel[2])
                fshift = fshift * (1 - n_d0t_v/C)
        
                doppler_shifted.append(fshift/1e3)
                time_shifted.append(dt.timedelta(seconds=rtime[final_n]) + t_time)
                
                allsec.append(rtime[final_n])
                alldop.append(fshift-freq)
                allthetas.append(thetas[ri])

                ray = r
            
            pulse_dop.append(doppler_shifted)
            pulse_tdelay.append(time_shifted)
        
        # last level
        tnt_dop.append(pulse_dop)
        tnt_t.append(pulse_tdelay)

        print('tnt time is', tim)

    return tnt_dop, tnt_t, alldop, allsec, allthetas, ray