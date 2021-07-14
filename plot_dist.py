import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from ray_plots import plotray2D
from raytracer_utils import read_rayfile
from constants_settings import *
from convert_coords import convert2
import datetime as dt
md = 6

ray_datenum = dt.datetime(2020,5,19,15,47,45,tzinfo=dt.timezone.utc)

ray_out_dir = '/media/rileyannereid/DEMETER/SR_output/2020-05-19_15_47_45/'
file_titles = os.listdir(ray_out_dir)
md_name = 'mode'+str(md)
raylist = []
for filename in file_titles:
    if '.ray' in filename and md_name in filename:
        if 'Main' in filename:
            pass
        elif '01' in filename or '02' in filename or '03' in filename:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))
            print(filename)

# plotting coords
crs = 'GEO'
carsph = 'car'
units = ['Re','Re','Re']

ray_coords = []
n_magsi = []
n_magsf = []
n_mags = []

thetas_save = []
nmax = plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', ['Re','Re','Re'], md, show_plot=False,plot_density=True,damping_vals=None)

for r in raylist:
    w = r['w']
    tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))

    # get nmag -- INITIAL AND all
    nmagi = np.sqrt(r['n'].x.iloc[0]**2 +r['n'].y.iloc[0]**2+r['n'].z.iloc[0]**2)
    nmag = list(np.sqrt(r['n'].x**2 +r['n'].y**2+r['n'].z**2))
    #if r['stopcond'] == 1:
    n_magsf.append(nmag[-1])
    n_magsi.append(nmagi)
    n_mags.append(nmag)

    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
    new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], 'GEO', 'sph', ['m','deg','deg'])
            
    # save it
    ray_coords.append(new_coords)
            
    B   =  r['B0'].iloc[0]
    Bmag = np.linalg.norm(B)
    bunit = B/Bmag

    # get the initial wna (already confirmed this is equiv to the initial)
    kveci = [(w/C)*r['n'].x.iloc[0],(w/C)*r['n'].y.iloc[0],(w/C)*r['n'].z.iloc[0]]
    kunit = np.array(kveci)/np.linalg.norm(kveci)
    alpha = np.arccos(np.dot(kunit, bunit))
    alphaedg = float(alpha)*R2D
    thetas_save.append(alphaedg)
#nmax=2045
#the_weight = (np.array(n_magsi)**2)/nmax
fig, ax = plt.subplots(1,1)
plt.style.use('seaborn-whitegrid')
ax.scatter(thetas_save,the_weight,marker='.', color='black')
ax.set_ylim([0,0.5])
plt.show()
plt.close()