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
from ray_plots import plotray2D, plotrefractivesurface, plot_damp, plotgeomfactor, stix_parameters

rayfile_directory = '/home/rileyannereid/workspace/SR_output' # store output here

# set mode and nrays and time increment
md = 7
time_inc = [0]
run_the_rays = False
nrays = 10000

# set burst start time and frequency
dd = dt.datetime(2020,4,6,22,5,0,tzinfo=dt.timezone.utc)
freq = 8.2e3

for tt in time_inc:
    ray_datenum = dd + dt.timedelta(seconds=tt)

    # we need the positions of the satellites -- use the sat class
    dsx = sat()             # define a satellite object
    dsx.catnmbr = 44344     # provide NORAD ID
    dsx.time = ray_datenum  # set time
    dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

    # propagate the orbit! setting sec=0 will give you just the position at that time
    dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

    vpm = sat()    
    vpm.catnmbr = 45120 
    vpm.time = ray_datenum - dt.timedelta(minutes=3) # back up to see full path
    vpm.getTLE_ephem()    

    vpm.propagatefromTLE(sec=60*6, orbit_dir='future', crs='GEO', carsph='sph', units=['Re','deg','deg'])

    # ray start position
    ray_start = dsx.pos

    # check direction of vpm
    thetas = []
    phis = []
    # vpm latitude
    checklat = vpm.pos[len(vpm.pos)//2][1]
    if checklat > 0:
        hemi_mult = 0
    else:
        hemi_mult = 180

    thetas.append(hemi_mult)
    for n in range(0,int(nrays)-1):
        thetas.append(0)
        phis.append(0)

    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
    th_file = ray_out_dir+'/thetas_'+str(md)+'.txt'

    if run_the_rays == True:
        # get directions and resonance cone
        directions, ra, thetas, phis = getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis, md, select_random=True)
        th_save = thetas

        # theta save has one extra entry for the reseonance cone
        th_save.append(ra*R2D)
        # save it for weighting
        np.savetxt(th_file,th_save)
    
        nrays = len(thetas) # how many rays -- THIS MUST BE EQUAL IN LENGTH TO THETAS AND PHIS
        
        positions = [ray_start[0] for n in range(nrays)]
        freqs = [freq for n in range(nrays)] 
        single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md)
    else:
        # damping already run
        file_titles = os.listdir(ray_out_dir)

        f=open(th_file)
        thetas_save = []
        for line in f:
            thetas_save.append(float(line))
        f.close()
        
        # extract ra
        # ra = thetas_save[-1]

        raylist = []
        for filename in file_titles:
            if '.ray' in filename and str(md) in filename:
                raylist += read_rayfile(os.path.join(ray_out_dir, filename))

        # get damping output to scale power
        """
        f = open(ray_out_dir+'/ray_out_mode'+str(md)+'_landau.txt')
        lines = []
        header = ['t','d']
        for line in f:
            if line.split() == header:
                continue
            else:
                lines.append(float(line.split()[1]))
        damping=lines
        
        for i in range(1893654 - 1829130):
            damping.append(1)
        """
        # plots relative power 2D, returns maximum refractive index from ra (uses Stix params function)
        nmax = plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md, show_plot=False,plot_density=True,damping_vals=None,theta_file=th_file)

        # plotting coords
        crs = 'GEO'
        carsph = 'car'
        units = ['Re','Re','Re']

        ray_coords = []
        n_magsi = []
        n_mags = []

        ray_count = 0
        # go through all rays, get the coords, n
        for r in raylist:
            ray_count += 1
            tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))

            # get nmag -- INITIAL AND all
            nmagi = np.sqrt(r['n'].x.iloc[0]**2 +r['n'].y.iloc[0]**2+r['n'].z.iloc[0]**2)
            nmag = list(np.sqrt(r['n'].x**2 +r['n'].y**2+r['n'].z**2))

            n_magsi.append(nmagi)
            n_mags.append(nmag)

            tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
            new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)
            
            # save it
            ray_coords.append(new_coords)

        # plot a histogram real quick of the rays that reflected
        """
        for rcs in ray_coords:
            ray_sph = convert2(rcs, tvec_datetime, crs, carsph, units, 'GEO', 'sph', ['m','deg','deg'])
            last_la = 0 
            last_lg = -1
            for rc in ray_sph:
                lg = np.abs(last_la) - np.abs(rc[1]) # latitude gradient 
                if np.sign(lg) != np.sign(last_lg):
                    print(rc[0]) # save the turn altitude
                last_lg = lg
                last_la = rc[1]
        """
        # set up plot
        fig = plt.figure(figsize=(11,8))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.coastlines()   

        # convert one ray at a time oops
        endlats = []
        endlons = []

        oopsrays = 0

        weights = []
        vlats = []
        vlons = []

        # this converts vpm longitude for plotting
        for v in vpm.pos:
            if v[2] > 180:
                long_new = v[2] - 180
                long_new = -1*(180 - long_new)
            else:
                long_new = v[2]
            vlons.append(long_new)
            vlats.append(v[1])
        checklon = vpm.pos[len(vpm.pos)//2][2]
        if checklon > 180:
            long_new = checklon - 180
            long_new = -1*(180 - long_new)
        else:
            long_new = checklon
        checklon = long_new

        # find where to center the plot
        extent = [checklon-5,checklon+5,checklat-5,checklat+5]
        damp_ind = 0

        # loop through all rays again
        for rr,nmagi,nmag in zip(ray_coords,n_magsi,n_mags):
            damp_ind += len(rr)-1

            # convert to spherical for plotting
            ray_sph = convert2(rr, tvec_datetime, crs, carsph, units, 'GEO', 'sph', ['m','deg','deg'])
            ra = [rs[0] for rs in ray_sph]
            rl = [rs[1] for rs in ray_sph]
            rlo = [rs[2] for rs in ray_sph]

            for ai, (alt_check, lat_check, lon_check) in enumerate(zip(ra,rl,rlo)):
                n_final = nmag[ai]

                if (alt_check-R_E) < 500e3 and np.sign(lat_check)==np.sign(checklat):
                    # fix longitude sign
                    if lon_check > 180:
                        long_new = lon_check - 180
                        long_new = -1*(180 - long_new)
                    else:
                        long_new = lon_check
                    endlats.append(lat_check)
                    endlons.append(long_new)
                    
                    # weight by initial wavenormal
                    the_weight = nmagi/nmax

                    # weight each ray to get relative power in each bin
                    # includes the final index of refraction as well
                    #weights.append(damping[damp_ind]*the_weight/n_final)
                    weights.append(the_weight/n_final)

                    break

        print('lost ', ray_count - len(endlats))
        print('rays found ',ray_count)

        binnum = 40
        binlon = np.linspace(extent[0],extent[1],num=binnum)
        binlat = np.linspace(extent[2],extent[3],num=binnum)

        # bin everything and find power density
        h = ax.hist2d(endlons,endlats,bins = [np.array(binlon),np.array(binlat)],weights=weights,alpha=0) #norm=mpl.colors.LogNorm(vmin=10e-2,vmax=nrays/10),zorder=10,alpha=0.5)

        # get area of each bin -- assume same area
        start_d = (binlat[0],binlon[0]) 
        end_d = (binlat[1],binlon[0])

        start_d2 = (binlat[1],binlon[0]) 
        end_d2 = (binlat[1],binlon[1])

        # find area
        a1 = geodesic(start_d,end_d).meters
        a2 = geodesic(start_d2,end_d2).meters
        bin_area = a1*a2

        # how many total watts? 
        wpr = 10/nrays
        # finally get to V/m 
        E_field = np.sqrt( 2*wpr*h[0] / ( C*EPS0*bin_area ))
        ef_c = ax.pcolormesh(np.array(binlon),np.array(binlat), E_field, zorder=10,alpha=0.5)

        cmin,cmax = h[3].get_clim()
        cax = ax.inset_axes([1.15, 0.2, 0.05, 0.6], transform=ax.transAxes)

        cbar = plt.colorbar(ef_c,ax=ax, cax=cax,shrink=0.45)

        # plot DSX fieldline footprint
        T = getBline(ray_start[0],ray_datenum,475)
        
        # repack
        T_repackx = T.x
        T_repacky = T.y
        T_repackz = T.z
        T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
        T_footpoint = T_repack[-1]
        
        T_convert = convert2([T_footpoint], [ray_datenum], 'SM','car',['Re','Re','Re'], crs, 'sph', ['Re','deg','deg'])
        ax.scatter(T_convert[0][2],T_convert[0][1],c='k',marker='*')

        # plot vpm
        ax.plot(vlons,vlats) 
        ax.scatter(checklon,checklat,c='r',marker='*')   
        ax.set_extent(extent)

        # clean up
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=2, color='gray', alpha=0, linestyle='--')
        plt.title(dt.datetime.strftime(ray_datenum,'%Y-%m-%d %H:%M:%S'))

        plt.savefig(ray_out_dir+'/'+dt.datetime.strftime(ray_datenum,'%Y%m%d_%H%M%S') + 'rayspotLEO_'+str(md)+'.png')
        plt.show()

        plt.close()