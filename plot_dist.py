import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 

from raytracer_utils import read_rayfile
from constants_settings import *
from convert_coords import convert2

th_file = '/home/rileyannereid/workspace/SR_output/2020-04-06 22:04:30/thetas_1.txt'

f=open(th_file)
thetas_save = []
for line in f:
    thetas_save.append(float(line))
f.close()

md = 1 

ray_datenum = dt.datetime(2020,4,6,22,4,30,tzinfo=dt.timezone.utc)

ray_out_dir = '/home/rileyannereid/workspace/SR_output/2020-04-06 22:04:30/other/'
file_titles = os.listdir(ray_out_dir)

raylist = []
for filename in file_titles:
    if '.ray' in filename and str(md) in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))
# plotting coords
crs = 'GEO'
carsph = 'car'
units = ['Re','Re','Re']

ray_coords = []
n_magsi = []
n_magsf = []

v_phs = []
for r in raylist:
    w = r['w']
    tmp_coords = list(zip(r['pos'].x, r['pos'].y, r['pos'].z))
    tmpv_coords = list(zip(r['vprel'].x, r['vprel'].y, r['vprel'].z))

    # get nmag -- INITIAL AND FINAL
    nmagi = np.sqrt(r['n'].x.iloc[0]**2 +r['n'].y.iloc[0]**2+r['n'].z.iloc[0]**2)
    nmagf = np.sqrt(r['n'].x.iloc[-1]**2 +r['n'].y.iloc[-1]**2+r['n'].z.iloc[-1]**2)
    n_magsi.append(nmagi)
    n_magsf.append(nmagf)

    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
    new_coords = convert2(tmp_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)
    new_vcoords = convert2(tmpv_coords, tvec_datetime, 'SM', 'car', ['m','m','m'], crs, carsph, units)
    
    # get final phase vel
    v_ph = np.sqrt(new_vcoords[-1][0]**2 + new_vcoords[-1][1]**2 + new_vcoords[-1][2]**2)

    # save it
    ray_coords.append(new_coords)
    v_phs.append(v_ph)

nmax = 2045
the_weight = np.array(n_magsi)/nmax
plt.style.use('seaborn-whitegrid')
plt.scatter(thetas_save[:-1],the_weight,marker='.', color='black')

plt.show()
plt.close()