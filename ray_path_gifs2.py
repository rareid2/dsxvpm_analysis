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
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap
# helpers
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter

def find_xy(p1, p2, z):

    x1, y1, z1 = p1
    x2, y2, z2 = p2
    if z2 < z1:
        return find_xy(p2, p1, z)

    # need to work on this
    x = np.interp(z, (z1, z2), (x1, x2))
    y = np.interp(z, (z1, z2), (y1, y2))

    return x, y

def rotate(points, normal, around):
  # Let's rotate the points such that the normal is the new Z axis

  old_x_axis = np.array([1, 0, 0])

  z_axis = normal
  y_axis = np.cross(old_x_axis, z_axis)
  x_axis = np.cross(z_axis, y_axis)
  
  axis = np.stack([x_axis, y_axis, z_axis])

  return np.dot(points - around, axis.T)

rayfile_directory = '/media/rileyannereid/DEMETER/SR_output' # store output here

# set burst start time and frequency
ray_datenum = dt.datetime(2020,6,1,1,46,50,tzinfo=dt.timezone.utc)
freq = 8.2e3
md = 7 # or 7
nrays=10000
checklat = 50 # BE SURE TO BE CHAGNING THIS
# lets chunk into a time vector
tnum =20
t = np.linspace(0,0.27,num=tnum)

"""
# this is all to get avg time --- comment this out
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
md = 7 # or 7
mode_name = 'mode' + str(md)
raylist = []
file_titles = os.listdir(ray_out_dir)
for filename in file_titles:
    if '.ray' in filename and mode_name in filename and '01' in filename:
        if 'Main' in filename:
            pass
        else:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))
            print(filename)
avg_time = 0
r_c=0
for r in raylist:
    if r['stopcond'] == 1:
        avg_time += max(r['time'])
        print(max(r['time']))
        r_c+=1
print('avg time', avg_time/r_c)

"""
# get the directory
ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')

# damping already run
file_titles = os.listdir(ray_out_dir)

mode_name = 'mode' + str(md)
raylist = []

# use mode name to avoid workers of the same label
for filename in file_titles:
    if '.ray' in filename and mode_name in filename:
        if 'Main' in filename:
            pass
        else:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))
            print(filename)

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
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter(T_repackx,T_repacky,T_repackz)
#ax.scatter(ray_start[0][0],ray_start[0][1],ray_start[0][2])
#plt.show()
#plt.close()

#for r in raylist:
#    print(r['stopcond'])
#    print(max(r['time']))

# loop through time
for si,s in enumerate(t):
    print('step', si, 'out of' , len(t))

    r_savex = []
    r_savey = []
    # go through all rays, get the coords, n
    ray_count = 0
    current_step = []
    mlats = []
    bad_ray = 0
    
    for ri,r in enumerate(raylist):
        ray_count += 1

        rt = np.array(r['time'])
        absolute_val_array = np.abs(rt - s)
        smallest_difference_index = absolute_val_array.argmin()

        closest_element = rt[smallest_difference_index]
        if np.abs(closest_element - s) > 0.1:
            #print(closest_element,s)
            continue 

        tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])

        # if altitude is increasing, ray is mirroring -- also confirm close (ish) to LEO
        mtime = None
        if r['stopcond'] != 1:
            # when did it mirror? 
            lasta = R_E
            loop_run = 0
            for ri, rc in enumerate(new_coords): 
                rc_sign = lasta - rc[0]
                # will be negative if altitide starts increasing..
                if np.sign(rc[1]) == np.sign(checklat) and np.sign(rc_sign) < 0 and ri > 1 and loop_run==0:
                    # 25km threshold? propagate to 
                    if rc[0]-R_E > 500e3:
                        #print((rc[0]-R_E)/1e3)
                        #when did they reflect? 
                        mtime = r['time'].iloc[ri]
                        loop_run = 1
                        continue
                    # go to the next ray
                # reset
                lasta = rc[0]
        
        if mtime:
            if mtime < s+0.01:
                bad_ray+=1
                continue # mirrored before this time step w some tolerance
            else:
                pass
                # ray mirrored but not yet 

        # get all the r's at that timestep -- if and only if the ray is going down!
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

        z = np.linspace(b_spot1[2],b_spot2[2],num=100)  
        new_pts = []
        for zz in z:
            new_pt = find_xy(b_spot1,b_spot2,zz)
            new_pt = [new_pt[0], new_pt[1], zz]
            new_pts.append(new_pt)
        b_zoom.append(new_pts)

        b_spot1 = T_repack[small_distances_idx[1]]
        b_spot2 = T_repack[small_distances_idx[2]]

        z = np.linspace(b_spot1[2],b_spot2[2],num=100)
        new_pts = []
        for zz in z:
            new_pt = find_xy(b_spot1,b_spot2,zz)
            new_pt = [new_pt[0], new_pt[1], zz]
            new_pts.append(new_pt)
        b_zoom.append(new_pts[1:])

        # unpack
        bline_final = [item for sublist in b_zoom for item in sublist]

        # quick check
        #fig = plt.figure()
        #ax = plt.axes(projection='3d')
        #ax.scatter(T_repackx,T_repacky,T_repackz)
        #bx = [bf[0] for bf in bline_final]
        #by = [bf[1] for bf in bline_final]
        #bz = [bf[2] for bf in bline_final]

        #ax.scatter(bx,by,bz)
        
        #plt.show()
        #plt.close()

        # NOW find the closest point and the index to the next one, point NORTH
        dist = [np.sqrt( (r_step[0] - b_pt[0])**2 + (r_step[1] - b_pt[1])**2 + (r_step[2] - b_pt[2])**2 ) for b_pt in bline_final]
        dist = np.array(dist)
        smallest_difference_index = dist.argmin()
        closest_element = bline_final[smallest_difference_index]
        nclosest_element = bline_final[smallest_difference_index-1]
        
        # grab mlat real quick
        cvals = coord.Coords([closest_element], 'SM', 'car')
        cvals.ticks = Ticktock([ray_datenum], 'ISO') # add ticks
        B_mag = cvals.convert('MAG', 'sph')
        mlat = B_mag.lati
        mlats.append(mlat)

        # finally, take normal
        Bunit = np.array(nclosest_element) - np.array(closest_element)
        Bunit = Bunit/np.linalg.norm(Bunit)

        # translate then rotate??
        #d = (Bunit[0]*closest_element[0] + Bunit[1]*closest_element[1] + Bunit[2]*closest_element[2])

        # another quick check -- north or south???
        #fig = plt.figure()
        #ax = plt.axes(projection='3d')
        #ax.scatter(T_repackx,T_repacky,T_repackz)
        #ax.scatter(closest_element[0],closest_element[1],closest_element[2],c='g')
        #ax.scatter(nclosest_element[0],nclosest_element[1],nclosest_element[2],c='b')
        #ax.quiver(closest_element[0],closest_element[1],closest_element[2],Bunit[0],Bunit[1],Bunit[2],length=10000)    
        #plt.show()
        #plt.close()

        r_step = np.array(r_step)

        out = rotate(r_step,Bunit,closest_element)

        # save all ending 'offsets'
        r_savex.append(out[0])
        r_savey.append(out[1])
    
    print('mirrored',bad_ray)

    #plt.scatter(r_savex,r_savey)
    #plt.show()
    #plt.close()

    fig, ax = plt.subplots(figsize=(10,12),tight_layout=True)
    
    # here is the bfield
    plt.scatter(0,0,marker="*",s=10,zorder=2)

    rx = np.array(r_savex)
    ry = np.array(r_savey)

    nan_array = np.isnan(rx)
    not_nan_array = ~ nan_array
    rxf = rx[not_nan_array]/1e3

    nan_array = np.isnan(ry)
    not_nan_array = ~ nan_array
    ryf = ry[not_nan_array]/1e3

    h = ax.hist2d(rxf,ryf,bins = 50,alpha=0)
    xedges = h[1]
    yedges = h[2]

    dx = xedges[1] - xedges[0]
    dy = yedges[1] - yedges[0]

    # set levels
    cbar_lvls = np.linspace(0,80,9)

    # Choose colormap which will be mixed with the alpha values
    cmap = pl.cm.coolwarm

    # Get the colormap colors
    my_cmap = cmap(np.arange(cmap.N))
    # Define the alphas in the range from 0 to 1
    alphas = 0.5*np.ones(cmap.N)
    # Define the background as white
    BG = np.asarray([1., 1., 1.,])
    # Mix the colors with the background
    for i in range(cmap.N):
        my_cmap[i,:-1] = my_cmap[i,:-1] * alphas[i] + BG * (1.-alphas[i])
    # Create new colormap which mimics the alpha values
    my_cmap = ListedColormap(my_cmap)

    ray_contour = ax.contourf(xedges[:-1]+dx/2,yedges[:-1] + dy/2, np.transpose(h[0]), cbar_lvls, zorder=1,cmap=my_cmap,extend='max')
    cbar = plt.colorbar(ray_contour,ax=ax,ticks=cbar_lvls,orientation='vertical',shrink=0.63,pad=0.015,label='# rays')
    
    plt.xlabel('km')
    plt.ylabel('km')

    left, right = plt.xlim()
    # which is bigger
    if np.abs(left) > np.abs(right):
        limi = np.abs(left)
    else:
        limi = np.abs(right)    
    
    plt.xlim([-limi,limi])
    plt.ylim([-limi,limi])

    ax.set_aspect('equal',adjustable='box')
    ax.set_facecolor('#a8b4ec')
    
    if limi < 1:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    elif limi < 1000:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # get average mlat
    m = 0
    for ml in mlats:
        m+=ml
    mlat_avg = m/len(mlats)
    mlat_str = str(int(mlat_avg))
    plt.title('mlat = ' + mlat_str+ ', # rays = '+ str(nrays-bad_ray))

    cname = ray_out_dir+ '/figures/mode'+ str(md)+'_raydist_'+str(si)+'.png'
    #plt.show()
    plt.savefig(cname,bbox_inches='tight')
    plt.close()
    imgs.append(cname)


simplegifs(imgs,ray_out_dir+ '/figures/'+ str(md)+'_raydist.gif')

# rotation
#denom = Bunit[0]**2 + Bunit[1]**2 + Bunit[2]**2
#costh = Bunit[2] / np.sqrt(denom)
#sinth = np.sqrt ( (Bunit[0]**2 + Bunit[1]**2) / denom )
#u1 = Bunit[1] / np.sqrt(Bunit[0]**2 + Bunit[1]**2)
#u2 = - Bunit[0] / np.sqrt(Bunit[0]**2 + Bunit[1]**2)

# rotation of plane normal to B rotated onto XY plane
#rot_m  = np.matrix([[costh+(u1**2)*(1-costh), u1*u2*(1-costh), u2*sinth],[u1*u2*(1-costh) , costh + (u2**2)*(1-costh), -u1*sinth],[-u2*sinth, u1*sinth, costh]])

# rotated B should now be at origin
# translate first
#B_translate = np.matrix([ [closest_element[0]],[closest_element[1]],[closest_element[2]-(d/Bunit[2])] ])
#r_step = np.matrix([ [r_step[0]],[r_step[1]],[r_step[2]-(d/Bunit[2])] ])

# finally, translate again
#r_final = np.matmul(rot_m, r_step)
#b_final = np.matmul(rot_m, B_translate)

#r_final_translate = [float(r_final[0]) - float(b_final[0]), float(r_final[1]) - float(b_final[1])]

# no
