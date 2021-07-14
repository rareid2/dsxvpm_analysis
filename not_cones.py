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
from planes import isect_line_plane_v3
from simplegifs_func import simplegifs
import math

# helpers
from scipy.interpolate import interp1d

def find_xy(p1, p2, z):

    x1, y1, z1 = p1
    x2, y2, z2 = p2
    if z2 < z1:
        return find_xy(p2, p1, z)

    # need to work on this
    x = np.interp(z, (z1, z2), (x1, x2))
    y = np.interp(z, (z1, z2), (y1, y2))

    return x, y

rayfile_directory = '/media/rileyannereid/DEMETER' # store output here

# set burst start time and frequency
dd = dt.datetime(2020,4,6,22,4,30,tzinfo=dt.timezone.utc)
freq = 8.2e3
md = 7 # or 7
ray_datenum = dd + dt.timedelta(seconds=30)

ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H_%M_%S')
th_file = ray_out_dir+'/thetas_'+str(md)+'.txt'

# damping already run
file_titles = os.listdir(ray_out_dir)

f=open(th_file)
thetas_save = []
for line in f:
    thetas_save.append(float(line))
f.close()

mode_name = 'mode' + str(md) + '.'
raylist = []
r_savex = []
r_savey = []

# use mode name to avoid workers of the same label
for filename in file_titles:
    if '.ray' in filename and mode_name in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))
        print(filename)

# lets chunk into a time vector
t = np.linspace(0,0.4,num=1)
t = [0.2]
imgs = []

# we need the positions of the satellites -- use the sat class
dsx = sat()             # define a satellite object
dsx.catnmbr = 44344     # provide NORAD ID
dsx.time = ray_datenum  # set time
dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

# propagate the orbit! setting sec=0 will give you just the position at that time
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

# now we have ray start point in the correct coordinates (SM cartesian in m)
ray_start = dsx.pos

# now get Bline!
T = getBline(ray_start[0], ray_datenum, 475)

# repack -- comes out in RE units
T_repackx = T.x * R_E
T_repacky = T.y * R_E
T_repackz = T.z * R_E
T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]

# quick check
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(T_repack)
ax.scatter(ray_start)
plt.show()
plt.close()

raylist = raylist[:10]
# loop through time
for si,s in enumerate(t):
    print(si)
        
    # go through all rays, get the coords, n
    ray_count = 0
    current_step = []
    
    for ri,r in enumerate(raylist):
        rt = np.array(r['time'])

        ray_count += 1

        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))

        # get all the r's at that timestep
        r_step = np.array(tmp_coords[smallest_difference_index])
        
        # now find the closest 3 points
        # distance to all points
        dist = [np.sqrt( (r_step[0] - b_pt[0])**2 + (r_step[1] - b_pt[1])**2 + (r_step[2] - b_pt[2])**2 ) for b_pt in T_repack]
        
        # closest 3
        dist_array = np.array(dist)

        # index of closest 3 points
        small_distances_idx = np.argsort(dist_array)[:3]
        small_distances_idx.sort()

        # build an array of 1000 on both sides
        b_zoom = []
        b_spot1 = T_repack[small_distances_idx[0]]
        b_spot2 = T_repack[small_distances_idx[1]]

        z = np.linspace(b_spot1[2],b_spot2[2],num=1000)  
        new_pts = []
        for zz in z:
            new_pt = find_xy(b_spot1,b_spot2,zz)
            new_pt = [new_pt[0], new_pt[1], zz]
            new_pts.append(new_pt)
        b_zoom.append(new_pts)

        b_spot1 = T_repack[small_distances_idx[1]]
        b_spot2 = T_repack[small_distances_idx[2]]

        z = np.linspace(b_spot1[2],b_spot2[2],num=1000)
        new_pts = []
        for zz in z:
            new_pt = find_xy(b_spot1,b_spot2,zz)
            new_pt = [new_pt[0], new_pt[1], zz]
            new_pts.append(new_pt)
        b_zoom.append(new_pts[1:])

        # unpack
        bline_final = [item for sublist in b_zoom for item in sublist]

        # quick check
        #plt.scatter(T_repackx,T_repackz)
        #new_ptx = [np[0] for np in bline_final]
        #new_ptz = [np[2] for np in bline_final]
        #plt.scatter(new_ptx,new_ptz)
        #plt.show()
        #plt.close()

        # NOW find the closest point and the index to the next one, point NORTH
        dist = [np.sqrt( (r_step[0] - b_pt[0])**2 + (r_step[1] - b_pt[1])**2 + (r_step[2] - b_pt[2])**2 ) for b_pt in bline_final]
        dist = np.array(dist)
        smallest_difference_index = dist.argmin()
        closest_element = bline_final[smallest_difference_index]
        # now we have closest element

        # lets just take a radial slice
        # how do we get XY distance tho? 

        nclosest_element = bline_final[smallest_difference_index-1]

        # finally, take normal
        Bunit = np.array(nclosest_element) - np.array(closest_element)
        Bunit = Bunit/np.linalg.norm(Bunit)

        # translate then rotate
        d = (Bunit[0]*closest_element[0] + Bunit[1]*closest_element[1] + Bunit[2]*closest_element[2])

        # another quick check -- north or south???
        #plt.scatter(T_repackx,T_repackz)
        #plt.scatter(closest_element[0],closest_element[2],c='g')
        #plt.scatter(nclosest_element[0],nclosest_element[2],c='b')
        #plt.quiver(closest_element[0],closest_element[2],Bunit[0],Bunit[2])    
        #plt.show()
        #plt.close()

        denom = Bunit[0]**2 + Bunit[1]**2 + Bunit[2]**2
        costh = Bunit[2] / np.sqrt(denom)
        sinth = np.sqrt ( (Bunit[0]**2 + Bunit[1]**2) / denom )
        u1 = Bunit[1] / np.sqrt(Bunit[0]**2 + Bunit[1]**2)
        u2 = - Bunit[0] / np.sqrt(Bunit[0]**2 + Bunit[1]**2)

        # rotation of plane normal to B rotated onto XY plane
        rot_m  = np.matrix([[costh+(u1**2)*(1-costh), u1*u2*(1-costh), u2*sinth],[u1*u2*(1-costh) , costh + (u2**2)*(1-costh), -u1*sinth],[-u2*sinth, u1*sinth, costh]])

        # rotated B should now be at origin
        # translate first
        B_translate = np.matrix([ [closest_element[0]],[closest_element[1]],[closest_element[2]-(d/Bunit[2])] ])
        r_step = np.matrix([ [r_step[0]],[r_step[1]],[r_step[2]-(d/Bunit[2])] ])

        # finally, translate again
        r_final = np.matmul(rot_m, r_step)
        b_final = np.matmul(rot_m, B_translate)

        r_final_translate = [float(r_final[0]) - float(b_final[0]), float(r_final[1]) - float(b_final[1])]

        # save all ending 'offsets'
        r_savex.append(r_final_translate[0])
        r_savey.append(r_final_translate[1])

    
    #plt.scatter(r_savex,r_savey)
    #plt.show()
    #plt.close()

    fig, ax = plt.subplots()

    rx = np.array(r_savex)
    ry = np.array(r_savey)

    nan_array = np.isnan(rx)
    not_nan_array = ~ nan_array
    rxf = rx[not_nan_array]

    nan_array = np.isnan(ry)
    not_nan_array = ~ nan_array
    ryf = ry[not_nan_array]

    h = ax.hist2d(rxf,ryf,bins = 50,alpha=0)
    xedges = h[1]
    yedges = h[2]

    dx = xedges[1] - xedges[0]
    dy = yedges[1] - yedges[0]
    
    ray_contour = ax.contourf(xedges[:-1]+dx/2,yedges[:-1] + dy/2, h[0], 10, zorder=10,alpha=0.5,cmap='coolwarm')
    cbar = plt.colorbar(ray_contour,ax=ax,orientation='horizontal',pad=0.03,label='# rays')
    plt.xlim([-50000,50000])
    plt.ylim([-50000,50000])

    #plt.clf()
    #plt.imshow(heatmap.T, extent=extent, origin='lower')
    #plt.xlim([-0.25,0.25])
    #plt.ylim([-0.25,0.25])
    cname = ray_out_dir+ '/'+ str(md)+'_raydist_'+str(si)+'.png'
    plt.savefig(cname)
    plt.close()
    imgs.append(cname)


imgs = [ray_out_dir+ '/'+ str(md)+'_raydist_'+str(si)+'.png' for si in range(1,20)]
simplegifs(imgs,ray_out_dir+ '/'+ str(md)+'_raydist.gif')
