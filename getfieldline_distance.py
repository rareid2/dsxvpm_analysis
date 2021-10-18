import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from bfield import getBline
from satellites import sat
import datetime as dt
from convert_coords import convert2
from constants_settings import *
import numpy as np 
from trace_fieldline import trace_fieldline_ODE_3D
import matplotlib.pyplot as plt 
import time 

"""
dates = [dt.datetime(2020,4,6,22,4,45, tzinfo=dt.timezone.utc),dt.datetime(2020,4,26,7,5,5, tzinfo=dt.timezone.utc),dt.datetime(2020,5,8,21,49,20, tzinfo=dt.timezone.utc),
dt.datetime(2020,5,10,2,38,31, tzinfo=dt.timezone.utc),dt.datetime(2020,5,16,12,4,15, tzinfo=dt.timezone.utc),dt.datetime(2020,5,19,15,47,45, tzinfo=dt.timezone.utc),
dt.datetime(2020,5,20,18,27,45, tzinfo=dt.timezone.utc),dt.datetime(2020,5,25,22,47,40, tzinfo=dt.timezone.utc),dt.datetime(2020,5,28,2,19,30, tzinfo=dt.timezone.utc),
dt.datetime(2020,5,29,22,43,5, tzinfo=dt.timezone.utc),dt.datetime(2020,6,1,1,46,50, tzinfo=dt.timezone.utc),dt.datetime(2020,6,3,13,46,25, tzinfo=dt.timezone.utc),
dt.datetime(2020,6,6,10,40,10, tzinfo=dt.timezone.utc),dt.datetime(2020,6,6,19,56,10, tzinfo=dt.timezone.utc),dt.datetime(2020,6,7,17,49,15, tzinfo=dt.timezone.utc),
dt.datetime(2020,6,16,13,0,30, tzinfo=dt.timezone.utc),dt.datetime(2020,6,17,15,37,20, tzinfo=dt.timezone.utc),dt.datetime(2020,6,18,22,53,45, tzinfo=dt.timezone.utc),
dt.datetime(2020,6,21,19,44,35, tzinfo=dt.timezone.utc),dt.datetime(2020,7,4,12,35,20, tzinfo=dt.timezone.utc),dt.datetime(2020,7,23,21,24,5, tzinfo=dt.timezone.utc),
dt.datetime(2020,7,25,0,4,50, tzinfo=dt.timezone.utc),dt.datetime(2020,7,27,20,53,55, tzinfo=dt.timezone.utc),dt.datetime(2020,8,8,15,9,5, tzinfo=dt.timezone.utc),
dt.datetime(2020,8,8,23,55,55, tzinfo=dt.timezone.utc),dt.datetime(2020,8,17,21,20,35, tzinfo=dt.timezone.utc),dt.datetime(2020,8,20,18,34,50, tzinfo=dt.timezone.utc),]
"""
# loop in time, at each time point find the distance
#dirs=[0,-1,0,-1,-1,0,0,-1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,-1,-1,0,-1,0,-1]
#burst_lens = [60,120,120,120,120,120,120,120,120,120,60,120,120,120,120,120,120,120,120,10,120,120,120,120,120,120,120]

dates = [dt.datetime(2020,6,1,1,46,50, tzinfo=dt.timezone.utc)]
dirs = [-1]
burst_lens = [80]
for ray_datenum1, dir, burst_len in zip(dates,dirs,burst_lens):

    transverse_distances = []
    fieldline_distances = []
    vpm_pos1 = []
    vpm_pos2 = []

    foot_pos1 = []
    foot_pos2 = []

    # set the start position
    dsx = sat()             # define a satellite object
    dsx.catnmbr = 44344     # provide NORAD ID
    dsx.time = ray_datenum1  # set time
    dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

    vpm = sat()    
    vpm.catnmbr = 45120 
    vpm.time = ray_datenum1  # back up to see full path
    vpm.getTLE_ephem()    

    vpm.propagatefromTLE(sec=burst_len-1, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

    # propagate the orbit! setting sec=0 will give you just the position at that time
    dsx.propagatefromTLE(sec=burst_len-1, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

    for i in range(0,burst_len,1):
        ray_datenum = ray_datenum1+dt.timedelta(seconds=i)
        #print(ray_datenum)

        # trace fieldline to exactly VPM's altitude
        vpm_GEO = convert2([vpm.pos[i]], [ray_datenum], 'SM','car',['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
        vpm_GEO_alt = vpm_GEO[0][0]/1000 # alt in km
        vpm_alt = vpm_GEO_alt-(R_E/1000)

        # get fieldline at this point
        T = getBline(dsx.pos[i],ray_datenum,vpm_alt)
        # repack
        T_repackx = T.x
        T_repacky = T.y
        T_repackz = T.z
        T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
        T_footpoint = T_repack[dir]

        T_convert = convert2([T_footpoint], [ray_datenum], 'SM','car',['Re','Re','Re'], 'GEO', 'sph', ['Re','deg','deg'])

        # find the start point
        lastdd = 0
        for ti,t_spot in enumerate(T_repack):
            next_spot = np.array(dsx.pos[i])/R_E # convert to RE
            dd = np.sqrt((next_spot[0]-t_spot[0])**2 + (next_spot[1]-t_spot[1])**2 + (next_spot[2]-t_spot[2])**2)
            if dd<lastdd:
                ti_save = ti
            lastdd = dd

        #print('min distance at', ti_save, np.array(dsx.pos[i])/R_E, T_repack[ti_save])
        
        distance = 0

        if dir == 0:
            for ti in range(ti_save,0,-1):
                if ti == len(T_repack)-1:
                    break
                else:
                    t_spot = T_repack[ti]
                    next_spot = T_repack[ti-1]
                    dd = np.sqrt((next_spot[0]-t_spot[0])**2 + (next_spot[1]-t_spot[1])**2 + (next_spot[2]-t_spot[2])**2)
                    distance +=dd
        else:
            for ti in range(ti_save,len(T_repack)+1):
                if ti == len(T_repack)-1:
                    break
                else:
                    t_spot = T_repack[ti]
                    next_spot = T_repack[ti+1]
                    dd = np.sqrt((next_spot[0]-t_spot[0])**2 + (next_spot[1]-t_spot[1])**2 + (next_spot[2]-t_spot[2])**2)
                    distance +=dd

        #print('fieldline distance of', distance*R_E/1000) # convert to km

        #plt.plot(T.x[ti_save:],T.z[ti_save:])
        #plt.scatter(dsx.pos[0][0]/R_E,dsx.pos[0][2]/R_E)
        #plt.show()
        #plt.close()

        # find transverse distance
        t_dist = np.sqrt(((vpm.pos[i][0]/R_E)-T_footpoint[0])**2+((vpm.pos[i][1]/R_E)-T_footpoint[1])**2+((vpm.pos[i][2]/R_E)-T_footpoint[2])**2)

        #print('transverse distance', t_dist*R_E/1000)
        #save
        transverse_distances.append(t_dist*R_E/1000)
        fieldline_distances.append(distance*R_E/1000)

        vpm_pos1.append(vpm_GEO[0][1])
        vpm_pos2.append(vpm_GEO[0][2])
        foot_pos1.append(T_convert[0][1])
        foot_pos2.append(T_convert[0][2])
    print(ray_datenum1)
    print('min fieldline distance is', min(np.array(fieldline_distances)), 'min transverse distance is', min(np.array(transverse_distances)))


    plt.plot(vpm_pos2,vpm_pos1,'b')
    plt.plot(foot_pos2,foot_pos1,'r')
    plt.show()
    plt.close()
    #time.sleep(30)