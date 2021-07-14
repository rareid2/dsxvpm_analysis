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

rayfile_directory = '/media/rileyannereid/DEMETER/SR_output' # store output here
# BE SURE TO CHECK THAT THE CORRECT STOP CONDS ARE IN THE RAYTRACER
# BE SURE TO CHANGE RAY PLOTS TO CORRECT WATTAGE AND HERE BELOW
# set mode and nrays and time increment
md = 6 # or 7
run_the_rays = False
run_damp = False

nrays = 10000
nworkers = 16

tt_increment = [0]

# set burst start time and frequency
if run_the_rays == True:
    for tti in tt_increment:
        ray_datenum = dt.datetime(2020,5,19,15,47,45,tzinfo=dt.timezone.utc)+dt.timedelta(seconds=tti)
        freq = 8.2e3

        # we need the positions of the satellites -- use the sat class
        dsx = sat()             # define a satellite object
        dsx.catnmbr = 44344     # provide NORAD ID
        dsx.time = ray_datenum  # set time
        dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

        # propagate the orbit! setting sec=0 will give you just the position at that time
        dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

        vpm = sat()    
        vpm.catnmbr = 45120 
        vpm.time = ray_datenum - dt.timedelta(minutes=5) # back up to see full path
        vpm.getTLE_ephem()    

        vpm.propagatefromTLE(sec=60*10, orbit_dir='future', crs='GEO', carsph='sph', units=['Re','deg','deg'])

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

        # create output directory
        ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
        th_file = ray_out_dir+'/thetas_'+str(md)+'.txt'

        # get directions and resonance cone
        directions, ra, thetas, phis = getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis, md, select_random=True)

        nrays = len(thetas) # how many rays
        print(nrays) # just verify this is right
        
        # break up the starting directions 
        directions_list = [directions[int(i * (nrays/nworkers)):int((i+1)*nrays/nworkers)] for i in range(nworkers)]
        print('RUNNING ', nrays, nrays/nworkers)

        # same freq and starting position for all
        freqs_list = [[freq for p in range(len(d))] for d in directions_list]
        positions = [[ray_start[0] for p in range(len(d))] for d in directions_list]

        tvec = [ray_datenum for n in range(nworkers)]
        directory_list = [rayfile_directory for i in range(len(tvec))]
        mds = [md for i in range(len(tvec))]

        #single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md)
        parallel_run_rays(tvec, positions, directions_list, freqs_list, directory_list, mds)

elif run_damp:
    ray_datenum = dt.datetime(2020,6,1,1,46,50,tzinfo=dt.timezone.utc) # just run damping once -- but going to need to add this in correctly down there
    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
    # get damping output to scale power -- run for each mode
    
    f = open(ray_out_dir+'/md'+str(md)+'_landau.txt')

    header = ['t','d']
    all_rays = []
    all_times = []
    new_ray = []
    new_times = []
    ray_count_damp = 0
    # this just parses the damp file
    for line in f:
        if line.split() == header:
            all_rays.append(new_ray)
            all_times.append(new_times)
            new_ray = []
            new_times = []
            ray_count_damp +=1
        else:
            new_ray.append(float(line.split()[1]))
            new_times.append(float(line.split()[0]))
    f.close()
    # going to try and get an average time vector to apply
    
    # find max time
    # strip first and last (should be empty) 
    all_rays = all_rays[1:-2]
    all_times = all_times[1:-2]
    max_t = max([tt[-1] for tt in all_times])

    # bad loop oops
    tt_array = np.linspace(0,max_t,50)
    avg_damp = []
    t_sum_last = 0
    t_count_last = 0
    # loop through each time
    for t in tt_array:
        t_sum = 0
        t_count = 0
        # loop all the rays
        for rt, rd in zip(all_times,all_rays):
            # loop all the timepoints inside the rays
            for rr, dd in zip(rt,rd):
                if np.abs(rr - t) < 1e-2:
                    if np.isnan(dd) == False:      
                        t_sum+=dd
                        t_count+=1
        if t_sum == 0 or t_count == 0: 
            t_sum = t_sum_last
            t_count = t_count_last
        avg_damp.append(t_sum/t_count)
        t_sum_last = t_sum
        t_count_last = t_count

    # final clean up for infs
    for adi,ad in enumerate(avg_damp):
        if ad < 1:
            pass
        else:
            avg_damp[adi] = avg_damp[adi-1]
    
    # now have an array or avg damping coeff at each time step for any ray
    a_file = open(ray_out_dir+'/md'+str(md)+'_landau_processed.txt', "w")
    for dtime,adamp in zip(tt_array,avg_damp):
        line = [str(dtime), str(adamp)]  # convert to strings
        a_file.write(' '.join(line) + "\n") 
        
    a_file.close()

# just plotting
else:
    for tti in tt_increment:
        ray_datenum = dt.datetime(2020,6,1,1,46,50,tzinfo=dt.timezone.utc)+dt.timedelta(seconds=tti)
        ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
        freq = 30e3

        vpm = sat()    
        vpm.catnmbr = 45120 
        vpm.time = ray_datenum - dt.timedelta(minutes=5) # back up to see full path
        vpm.getTLE_ephem()    

        vpm.propagatefromTLE(sec=60*10, orbit_dir='future', crs='GEO', carsph='sph', units=['Re','deg','deg'])

        # we need the positions of the satellites -- use the sat class
        dsx = sat()             # define a satellite object
        dsx.catnmbr = 44344     # provide NORAD ID
        dsx.time = ray_datenum  # set time
        dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

        # propagate the orbit! setting sec=0 will give you just the position at that time
        dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])
        ray_start = dsx.pos

        # vpm latitude
        checklat = vpm.pos[len(vpm.pos)//2][1]
        if checklat > 0:
            hemi_mult = 0
        else:
            hemi_mult = 180

        # damping already run
        file_titles = os.listdir(ray_out_dir)

        raylist = []
        mode_name = 'mode' + str(md)
        # use mode name to avoid workers of the same label
        x = sorted(file_titles)
        for filename in x:
            if '.ray' in filename and mode_name in filename:
                if 'Main' in filename: # avoid call from bfield
                    pass
                else:
                    raylist += read_rayfile(os.path.join(ray_out_dir, filename))
                    print(filename)

        # get damping output to scale power -- run for each mode
        f = open(ray_out_dir+'/md'+str(md)+'_landau_processed.txt')
        tt_array = []
        avg_damp = []
        for line in f:
            lines = line.split()
            tt_array.append(float(lines[0]))
            avg_damp.append(float(lines[1]))
        f.close()

        # for plot ray 2D
        damping_list = [list(tt_array),avg_damp]

        # plots relative power 2D, returns maximum refractive index from ra (uses Stix params function)
        nmax = plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md, show_plot=False,plot_density=True,damping_vals=damping_list)

        ray_coords = []
        n_magsi = []
        n_mags = []

        ray_count = 0
        # go through all rays, get the coords, n
        thetas_save = []
        for r in raylist:
            ray_count += 1
            tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))

            # get nmag -- INITIAL AND all
            nmagi = np.sqrt(r['n'].x.iloc[0]**2 +r['n'].y.iloc[0]**2+r['n'].z.iloc[0]**2)
            nmag = list(np.sqrt(r['n'].x**2 +r['n'].y**2+r['n'].z**2))

            n_magsi.append(nmagi)
            n_mags.append(nmag)

            tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
            new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
            
            # save it
            ray_coords.append(new_coords)
            
            w = r['w']
            B   =  r['B0'].iloc[0]
            Bmag = np.linalg.norm(B)
            bunit = B/Bmag

            # get the initial wna (already confirmed this is equiv to the initial)
            kveci = [(w/C)*r['n'].x.iloc[0],(w/C)*r['n'].y.iloc[0],(w/C)*r['n'].z.iloc[0]]
            kunit = np.array(kveci)/np.linalg.norm(kveci)
            alpha = np.arccos(np.dot(kunit, bunit))
            alphaedg = float(alpha)*R2D
            thetas_save.append(alphaedg)

        # PLOT SET UP
        fig = plt.figure(figsize=(14,7),constrained_layout=True)
        grid = fig.add_gridspec(2, 3, height_ratios=[0.2,1], figure=fig)

        ax1 = fig.add_subplot(grid[0, :])
        ax2 = fig.add_subplot(grid[1, 0], projection=ccrs.PlateCarree())
        ax3 = fig.add_subplot(grid[1, 1], projection=ccrs.PlateCarree())
        ax4 = fig.add_subplot(grid[1, 2], projection=ccrs.PlateCarree())

        # plot clean up
        plt.rcParams['axes.edgecolor']='#333F4B'
        #plt.rcParams['axes.linewidth']=0.8
        plt.rcParams['xtick.color']='#333F4B'
        plt.rcParams['ytick.color']='#333F4B'
        
        # for the top row
        ax2.coastlines() 
        ax3.coastlines()
        ax4.coastlines()

        # convert one ray at a time oops
        endlats = []
        endlons = []

        oopsrays = 0

        weights = []
        wna_weights = []

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
        spatial_extent = 15
        extent = [checklon-spatial_extent,checklon+spatial_extent,checklat-spatial_extent,checklat+spatial_extent]

        # loop through all rays again lol
        for ray,rr,nmagi,nmag,th in zip(raylist,ray_coords,n_magsi,n_mags,thetas_save):
            ray_damp = 0

            ra = [rs[0] for rs in rr]
            rl = [rs[1] for rs in rr]
            rlo = [rs[2] for rs in rr]

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

                    wna_weights.append(th)
                    
                    # weight by initial wavenormal
                    the_weight = nmagi/nmax

                    # weight each ray to get relative power in each bin
                    # includes the final index of refraction as well
                    tt = list(ray['time'])
                    end_time = tt[-1]
                    for ti,ta in enumerate(tt_array):
                        if np.abs(end_time-ta) < 1e-1: # find which time point is closest (ish)
                            ray_damp = avg_damp[ti]

                    weights.append(ray_damp*the_weight/n_final)

                    break

        print('lost ', ray_count - len(endlats))
        print('rays found ',ray_count)

        binnum = 40
        binlon = np.linspace(extent[0],extent[1],num=binnum)
        binlat = np.linspace(extent[2],extent[3],num=binnum)

        # bin everything and find power density
        # hw has weights for power
        # nrays, wnas, power
        h = ax2.hist2d(endlons,endlats,bins = [np.array(binlon),np.array(binlat)],alpha=0) #norm=mpl.colors.LogNorm(vmin=10e-2,vmax=nrays/10),zorder=10,alpha=0.5)
        hh = ax3.hist2d(endlons,endlats,bins = [np.array(binlon),np.array(binlat)],weights=wna_weights,alpha=0) #norm=mpl.colors.LogNorm(vmin=10e-2,vmax=nrays/10),zorder=10,alpha=0.5)
        hw = ax4.hist2d(endlons,endlats,bins = [np.array(binlon),np.array(binlat)],weights=weights,alpha=0) #norm=mpl.colors.LogNorm(vmin=10e-2,vmax=nrays/10),zorder=10,alpha=0.5)

        # get area of each bin -- assume same area
        start_d = (binlat[0],binlon[0]) 
        end_d = (binlat[1],binlon[0])

        start_d2 = (binlat[1],binlon[0]) 
        end_d2 = (binlat[1],binlon[1])

        # find area
        a1 = geodesic(start_d,end_d).meters
        a2 = geodesic(start_d2,end_d2).meters
        bin_area = a1*a2

        dy = binlat[1] - binlat[0]
        dx = binlon[1] - binlon[0]

        # how many total watts? 
        wpr = 10/nrays
        # finally get to V/m 
        E_field = np.sqrt( 2*wpr*hw[0] / ( C*EPS0*bin_area ))
        # for power
        #ef_power = ax4.pcolormesh(np.array(binlon),np.array(binlat), E_field*10**6, zorder=10,alpha=0.5,cmap='coolwarm')

        # wavenormals -- old
        #bin_dens = np.zeros_like(h[0])
        #wna_dens = np.zeros_like(h[0])
        #for eli,(elo,ela) in enumerate(zip(endlons,endlats)):
        #    for ii, (bl1, bl2) in enumerate(zip(binlon, binlon[1:])):
        #        for kk, (al1, al2) in enumerate(zip(binlat, binlat[1:])):
        #            if bl1 < elo < bl2 and al1 < ela < al2:
        #                bin_dens[ii,kk] += 1
        #                wna_dens[ii,kk] += wna_weights[eli]

        #wna_avg = np.divide(wna_dens, bin_dens, out=np.zeros_like(wna_dens), where=bin_dens!=0)
        wna_avg = np.divide(hh[0], h[0], out=np.zeros_like(hh[0]), where=h[0]!=0)
        #with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
            #m3 = ax3.pcolormesh(ybins, xbins, sums / counts, cmap='coolwarm')
        #ef_ray = ax2.pcolormesh(np.array(binlon),np.array(binlat), h[0], zorder=10,alpha=0.5,cmap='coolwarm')
        #ef_wna = ax3.pcolormesh(np.array(binlon),np.array(binlat), wna_avg, zorder=10,alpha=0.5,cmap='coolwarm')
        
        # change to center points for contour
        binlat = np.array(binlat)
        binlon = np.array(binlon)

        ef_ray = ax2.contourf(binlon[:-1]+dx/2,binlat[:-1] + dy/2, np.transpose(h[0]), 10, zorder=10,alpha=0.5,cmap='coolwarm')
        ef_wna = ax3.contourf(binlon[:-1]+dx/2,binlat[:-1] + dy/2, np.transpose(wna_avg), 10, zorder=10,alpha=0.5,cmap='coolwarm')
        ef_power = ax4.contourf(binlon[:-1]+dx/2,binlat[:-1] + dy/2, np.transpose(E_field*10**6),10, zorder=10,alpha=0.5,cmap='coolwarm')

        cbar = plt.colorbar(ef_ray,ax=ax2,orientation='horizontal',pad=0.03,label='# rays')
        cbar = plt.colorbar(ef_wna,ax=ax3,orientation='horizontal',pad=0.03,label='deg')
        cbar = plt.colorbar(ef_power,ax=ax4,orientation='horizontal',pad=0.03,label='uV/m')

        # plot DSX fieldline footprint
        T = getBline(ray_start[0],ray_datenum,475)
        
        # repack
        T_repackx = T.x
        T_repacky = T.y
        T_repackz = T.z
        T_repack = [[tx, ty, tz] for tx,ty,tz in zip(T_repackx, T_repacky, T_repackz)]
        T_footpoint = T_repack[0]
        
        T_convert = convert2([T_footpoint], [ray_datenum], 'SM','car',['Re','Re','Re'], 'GEO', 'sph', ['Re','deg','deg'])
        
        # final plot stuff
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        for ax in [ax2,ax3,ax4]:
            ax.plot(vlons,vlats) 
            ax.scatter(checklon,checklat,c='r',marker='*')   
            ax.set_extent(extent)
            ax.scatter(T_convert[0][2],T_convert[0][1],c='k',marker='*')
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=2, color='gray', alpha=0, linestyle='--')
            gl.top_labels = True
            if ax == ax2:
                gl.left_labels = True
            else:
                gl.left_labels = False
            gl.right_labels =False
            gl.bottom_labels = False

            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER

        # lastly, the top histogram
        reflected_above = 0
        xaxis = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85]
        refl_th = np.zeros(len(xaxis))
        for rci,(rcs,th) in enumerate(zip(ray_coords,thetas_save)):
            lasta = R_E
            loop_run = 0
            for ri, rc in enumerate(rcs): 
                rc_sign = lasta - rc[0]
                # will be negative if altitide starts increasing..
                if np.sign(rc[1]) == np.sign(checklat) and np.sign(rc_sign) < 0 and ri > 1 and loop_run==0:
                    # 25km threshold? propagate to 
                    if rc[0]-R_E > 500e3:
                        #print((rc[0]-R_E)/1e3)
                        reflected_above+=1
                        # when did they reflect? 
                        #print(raylist[rci]['time'].iloc[ri])

                        # sort by initial WNA
                        for xi,xbar in enumerate(xaxis):
                            if np.abs(th) < xbar:
                                refl_th[xi] += 1
                                break
                        loop_run = 1
                        continue
                    # go to the next ray
                # reset
                lasta = rc[0]

        bar_heights = list(100*refl_th/reflected_above)
        print(reflected_above)

        ax1.vlines(x=xaxis, ymin=0, ymax=bar_heights, color='#9ca4dc', alpha=0.2, linewidth=5)
        ax1.plot(xaxis, bar_heights, "o", markersize=5, color='#9ca4dc', alpha=0.6)
        ax1.set_ylabel('% of ' + str(round(100*reflected_above/nrays,1))+'% \n that mirrored', fontsize=8, fontweight='black', color = '#333F4B')
        ax1.annotate('wna', (0.97, -0.16), xycoords='axes fraction', va='center')

        grid.update(wspace=0.01, hspace=0.05) # set the spacing between axes. 

        # finally title and save    
        if md == 6:
            strmd = 'GCPM'
        elif md == 7:
            strmd = 'Diff Eq'

        fig.suptitle(dt.datetime.strftime(ray_datenum,'%m %d %H:%M:%S') + ' in '+ strmd +' at ' + str(round(freq/1e3,1)) + ' kHz')
        plt.savefig(ray_out_dir+'/figures/'+dt.datetime.strftime(ray_datenum,'%Y%m%d_%H%M%S') + 'rayspotLEO_'+str(md)+'.png',bbox_inches='tight')
        #plt.show()

        plt.close()