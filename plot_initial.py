import numpy as np
import datetime as dt
import os
os.chdir('/home/rileyannereid/workspace/SR_interface')
import sys
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from constants_settings import *
from convert_coords import convert2
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
import PyGeopack as gp
import random
import matplotlib.pyplot as plt
import math
from planes import isect_line_plane_v3
rayfile_directory = '/home/rileyannereid/workspace/SR_output' # store output here

# set burst start time and frequency
dd = dt.datetime(2020,4,6,22,4,30,tzinfo=dt.timezone.utc)
freq = 8.2e3
md = 7 # or 7
ray_datenum = dd + dt.timedelta(seconds=30)

ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
th_file = ray_out_dir+'/thetas_'+str(md)+'.txt'
    
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []
for filename in file_titles:
    if '.ray' in filename and str(md) in filename and str('Main') in filename:
        print(filename)
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

nrays = 10000

# get b direction for this ray
for r in raylist:
    B0 = [r['B0'].x[0], r['B0'].y[0], r['B0'].z[0]]
    # vec in T in SM car coordinates
    # create unit vector 
    Bmag = np.sqrt(r['B0'].x[0]**2 + r['B0'].y[0]**2 + r['B0'].z[0]**2)
    Bunit = [r['B0'].x[0]/Bmag, r['B0'].y[0]/Bmag, r['B0'].z[0]/Bmag]

# also return resonance angle, can be useful for initializing rays
from ray_plots import stix_parameters
R, L, P, S, D = stix_parameters(r, 0, r['w']) # get stix params for initial time point
resangle = np.arctan(np.sqrt(-P/S))

converted_dirs = []
hemi_mult = 0
thetas = []
phis = []
resangle_deg = resangle *180/np.pi

for n in range(0,nrays):
    # sample theta as concentric circles around the z axis, max at resonance angle
    thetas.append((random.random()*(resangle_deg-3)))
    # uniform azimuth around the z axis
    phis.append(random.random()*360)

if Bunit[0] == 0 or Bunit[2] == 0:
    r1 = [1,(-1*Bunit[0]-Bunit[2])/Bunit[1],1]
else:
    r1 = [1,1,(-1*Bunit[1]-Bunit[0])/Bunit[2]]

r1 = np.array(r1)/np.linalg.norm(np.array(r1))
r2  = np.cross(r1,Bunit)
T_rotate = np.column_stack((r1,r2,Bunit))


for th,ph in zip(thetas,phis):
    r = 1/(np.cos(th*D2R))
    cone_vec = np.array([r*np.sin(th*D2R)*np.cos(ph*D2R),r*np.sin(th*D2R)*np.sin(ph*D2R),r*np.cos(th*D2R)]) 
    cone_vec = np.matmul(T_rotate,np.transpose(cone_vec))
    if hemi_mult == 180:
        zsign = -1
    else:
        zsign = 1

    cone_vec = cone_vec/np.linalg.norm(cone_vec)
    converted_dirs.append(zsign*cone_vec)

# find the plane
denom = Bunit[0]**2 + Bunit[1]**2 + Bunit[2]**2
costh = Bunit[2] / np.sqrt(denom)
sinth = np.sqrt ( (Bunit[0]**2 + Bunit[1]**2) / denom )
u1 = Bunit[1] / np.sqrt(denom)
u2 = - Bunit[2] / np.sqrt(denom)
rot_m  = np.array([[costh+u1**2*(1-costh), u1*u2*(1-costh), u2*sinth],[u1*u2*(1-costh) , costh + u2**2*(1-costh), -u1*sinth],[-u2*sinth, u1*sinth, costh]])

npp = np.matmul(rot_m,np.array(Bunit))
bx = npp[0]
by = npp[1]

pts = []
for cd in converted_dirs:

    if math.isnan(cd[0]):
        pass
    else:
        pt = isect_line_plane_v3([0,0,0],cd,Bunit,Bunit)
        pp = np.array(pt)
        npp = np.matmul(rot_m,pp)
        pts.append(npp)

# unpack and translate
px = [p[0]-bx for p in pts]
py = [p[1]-by for p in pts]

heatmap, xedges, yedges = np.histogram2d(px, py, bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower',zorder=100)
print(np.cos(resangle))
plt.Circle(( 0 , 0 ), np.cos(resangle),fill=False,zorder=101)
plt.colorbar()
#plt.xlim([-0.25,0.25])
#plt.ylim([-0.25,0.25])
cname = ray_out_dir+ '/'+ str(md)+'_raydist_initial.png'
plt.savefig(cname)
plt.close()