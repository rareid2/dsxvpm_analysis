# AT64ThCh.py
#
# Implementation of an diffusive-equilibrium density model for
# electrons, oxygen and hydrogen ions, based on Angerami and Thomas
# [1964] (AT64), Thomson [1987] and Chiu [1979].  The basic treatment
# follows AT64, with the definition of a "temperature modified
# geopotential height (z)".  Thomson includes an oxygen frictional
# term which results in a Chapman-like ionosphere.  Thomson's form of
# the temperature profile is also used here.  Chiu's inhomegeneous
# magnetic force term is also included.
#
# Carpenter-Anderson plasmasphere model (equatorial profile), with
# modifications due to Sheeley (plasmatrough), Moldwin (plasmapause
# location), Gallager (pp width) and Denton (latitudinal dependence).
#
# Default input parameters are chosen to try to emulate the
# Chiu/Abel&Thorne result.
#
# Assumes an eccentric dipole field model

# Michael Starks     AFRL/RVBX      16 Feb 2021


"""A diffusive equilibrium density model stitched to a modified Carpenter-Anderson plasmasphere model"""

# =======================
# Import required modules
# =======================
import numpy as np
import sys
import os
import datetime as dt 
import matplotlib.pyplot as plt
os.chdir('/home/rileyannereid/workspace/SR_interface')
sys.path.append('/home/rileyannereid/workspace/SR_interface') 
from bfield import getBdir, getBline
from convert_coords import convert2


def get_edens(rgeo, ray_datenum): # input in Mm
    M_P = 1.6726219e-27

    R_E = 6.3712                # Mm

    Omass = 2.67616e-26 # kg
    Hmass = 6.6904e-27 # kg

    # defaults
    PeakHeight=300
    ETransitionDensity=2.0e11
    OHTransitionHeight=400
    OHTransitionTemp=750
    NeutralTemp=1000
    TempGradient=800
    Kp=1

    peakheight = PeakHeight / 1.0e3
    etransdens = ETransitionDensity
    ohtransheight = OHTransitionHeight / 1.0e3
    ohtranstemp = OHTransitionTemp
    neutraltemp = NeutralTemp
    tempgradient = TempGradient
    kp = Kp

    # Calculate the plasmapause location (CA) and width (Gallagher nominal/dawn-dusk)
    Lpp = 5.6 - 0.46 * kp
    Lw = 0.14

    """Returns the electron and individual ion densities at a geocentric Cartesian location [Mm]"""

    # Initialize some of the outputs
    edens = 1.0         # to facilitate log-scale plotting (vice zero)
    idens = np.array((0.5,0.5), ndmin=2)
    imass = np.array((Omass, Hmass), ndmin=2)

    # Geocentric radial distance [Mm].  Abort if below 30 km
    r = np.linalg.norm(rgeo)

    # variation of gravity -NO
    g = 9.80665  # m/s^2
    #g = g0*(R_E/r) # m/s^2

    mpg = M_P * g      # proton mass times gravitational acceleration [kg m/(s^2)]

    if (r < (R_E + 0.030)):
        edens = 0.0
        idens = np.array((6, 6), ndmin=2)
        edens = 6
        return edens, idens

    phi = np.arcsin(rgeo[3]/r)
    L = (r/R_E) / (np.cos(phi)**2)

    h = r - R_E                         # height [Mm]

    r0 = R_E + ohtransheight           # geocentric transition height [Mm]
    
    R = r / r0                         # geocentric height ratio

    # Compute zbrat, the inhomogeneous field component of the z (Chiu
    # contribution); the magnetic field strength at the base level
    # (transition height) must be computed on *this* field-line.  The
    # computation here implicitly assumes a dipole field, but doesn't
    # make any assumptions about tilt or offset.  We use an eccentric
    # dipole.  When pulled through the equations, this becomes a
    # multiplicative factor

    # to find location at start alt, need to extrapolate
    start_alt = (ohtransheight+R_E)/(R_E) # now in Re

    if start_alt < (r/R_E):
        
        tvec_datetime = [ray_datenum]

        newt = convert2([rgeo*1e6], tvec_datetime, 'GEO', 'car', ['m','m','m'], 'SM', 'car', ['m','m','m'])
        T = getBline(newt[0], tvec_datetime[0], 400) # 400 km is the transition height
        
        # repack
        B_repackx = T.Bx
        B_repacky = T.By
        B_repackz = T.Bz
        B_repack = [[tx, ty, tz] for tx,ty,tz in zip(B_repackx, B_repacky, B_repackz)]

        B0 = np.linalg.norm(B_repack[-1])

        import PyGeopack as gp
        
        bdate = (ray_datenum.year * 1e4) + (ray_datenum.month*1e2) + ray_datenum.day
        ut = ray_datenum.hour + ray_datenum.minute/60 + ray_datenum.second/3600

        rgeo_bpos = np.array(newt)/(R_E*1e6)
        rgeo_bpos = rgeo_bpos[0]
        Bx,By,Bz = gp.ModelField(rgeo_bpos[0],rgeo_bpos[1],rgeo_bpos[2],bdate,ut,Model='T96',CoordIn='GSM',CoordOut='GSM')
        B = np.linalg.norm([Bx,By,Bz])

        zbrat = B/B0

    else:
        zbrat = 1
    """
    Babs, _Bdir = EccDipole.LocalField(rgeo)
    rmag = EccDipole.GEOtoMag(rgeo)
    rmag_sph = Cart2Sph(rmag)
    L = EccDipole.ComputeL(rgeo)
    s2 = np.square(np.cos(np.radians(rmag_sph[0,1])))
    arg = np.min((1.0, np.sqrt(s2 / R)))
    T0 = np.arcsin(arg)
    rmag0_sph = np.empty((1,3))
    rmag0_sph[0,2] = rmag_sph[0,2]
    rmag0_sph[0,1] = 90.0 - np.degrees(T0)
    rmag0_sph[0,0] = r0
    rmag0 = Sph2Cart(rmag0_sph)
    rgeo0 = EccDipole.MagtoGEO(rmag0)
    Babs0, _Bdir0 = EccDipole.LocalField(rgeo0)
    zbrat = Babs / Babs0
    """
    print(zbrat)
    # Compute the diffusive equilibrium elements
    a = tempgradient * r0 / ohtranstemp - 1.0                   # constant with height
    
    tt = (R * (1.0 + a) - a) / R
    
    zg = r0 / a * np.log(tt)                                                # Gravitational term of z
    Rp = (R_E + peakheight) / r0 	                # Geocentric height ratio of F2 peak

    c = 1.0 / ((Rp * (1.0 + a) - a) * Rp)	                                # constant with height
    H0 = 1.380658e-23 * neutraltemp / (16.0 * mpg) / 1.0e6      # Neutral O scale height [Mm]
    z = zg + c * H0 * np.exp((peakheight - h) / H0)                   # Modified geopotential height [Mm]
    T = ohtranstemp * tt                                              # Ion (electron) temperature [K]
    H1 = 1.380658e-23 * ohtranstemp / mpg / 1.0e6               # Hydrogen ion scale height [Mm]
    H3 = 1.380658e-23 * ohtranstemp / (16.0 * mpg) / 1.0e6      # Oxygen ion scale height [Mm]
    n10 = 0.5 * etransdens                                            # H+ concentration at transition height [m-3]
    n30 = 0.5 * etransdens	                                        # O+ concentration at transition height [m-3]

    # The predicted electron density
    ne = np.sqrt((etransdens * ohtranstemp) * zbrat * ((n10 * ohtranstemp) * np.exp(-z / H1)  + (n30 * ohtranstemp) * np.exp(-z / H3))) / T

    # Electron-ion density ratio
    R13 = (n10 / n30) * np.exp(z * ((H1 - H3) / (H1 * H3)))
    
    # Compute plasmatrough (Sheeley mean over local-time)
    SN = 124.0 * np.power((3.0 / L), 4.0) * 1.0e6

    # Plasmapause transition
    tran = 0.5 * np.tanh(3.4534 * (L - Lpp) / Lw ) + 0.5
    ne = (1.0 - tran) * ne  +  tran * SN

    # Complete outputs
    edens = ne
    idens = np.empty((1,2))

    idens[0,0] = ne / (1.0 + R13)
    if (R13 > 1e-30):
        idens[0,1] = ne / (1.0 + 1.0 / R13)
    else:
        idens[0,1] = 0.0
    #print(idens)
    return np.log10(edens), np.log10(idens)


# ----------------------- some test runs ---------------------
# start 30km out
R_E = 6.3712                # Mm

ray_datenum = dt.datetime(2020,8,25,9,36, tzinfo=dt.timezone.utc)
"""
xrego = np.linspace(-4e4/1e3,4e4/1e3,num=100)
zrego = np.linspace(-4e4/1e3,4e4/1e3,num=100)

edens = np.zeros((len(xrego),len(zrego)))
for i,x in enumerate(xrego):
    for j,z in enumerate(zrego):
        rgeo = np.array([x,0,z]) 
        edens[i,j], idens = get_edens(rgeo,ray_datenum)
    print(i)

fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
from matplotlib import cm
cmap = cm.get_cmap('viridis', 256*2)
psm = ax.pcolormesh(xrego,zrego,np.transpose(edens), cmap='turbo', vmin=6,vmax=13.75, rasterized=True)
fig.colorbar(psm, ax=ax)
ax.set_aspect('equal', 'box')
plt.show()
plt.close()
"""
"""
xrego = np.linspace(R_E,4e4/1e3,num=100)
magcrs = np.linspace(R_E*1e6,7*R_E*1e6,num=100)
magcrs = xrego
edens = np.zeros(len(magcrs))
odens = np.zeros(len(magcrs))
hdens = np.zeros(len(magcrs))

for i,x in enumerate(magcrs):
    maggeo = [[x,0,0]]
    rgeo = convert2(maggeo,[ray_datenum],'MAG','sph',['m','deg','deg'],'GEO','car',['m','m','m'])
    rgeo = [rgeo[0][0]/1e6,rgeo[0][1]/1e6,rgeo[0][2]/1e6]
    
    rgeo = maggeo[0]
    edens[i], idens = get_edens(np.array(rgeo),ray_datenum)
    odens[i] = idens[0][0]
    hdens[i] = idens[0][1]

plt.plot(xrego/R_E,edens,'b',label='e -')
plt.plot(xrego/R_E,odens,'r',label='O+')
plt.plot(xrego/R_E,hdens,'g',label='p+')
plt.legend()
plt.xlabel('Earth radii')
plt.ylabel('log(Ne [e/m^3])')
plt.ylim([6,12])
plt.grid()
plt.show()
#plt.xlabel()
plt.close()

#Lshell 3
newt = convert2([[3*R_E*1e6,0,0]], [ray_datenum], 'GEO', 'car', ['m','m','m'], 'SM', 'car', ['m','m','m'])
T = getBline(newt[0], ray_datenum, 400) # 400 km is the transition height

# repack
B_repackx = T.x
B_repacky = T.y
B_repackz = T.z
B_repack = [[tx, ty, tz] for tx,ty,tz in zip(B_repackx, B_repacky, B_repackz)]
tvec_datetime = [ray_datenum for tx in B_repack]
crs_l = convert2(B_repack, tvec_datetime,'SM','car',['Re','Re','Re'],'GEO','car',['m','m','m'])

magcrs = crs_l
edens = np.zeros(len(magcrs))
odens = np.zeros(len(magcrs))
hdens = np.zeros(len(magcrs))

for ci,cr in enumerate(crs_l):
    #print(cr)
    edens[ci], idens = get_edens(np.array(cr)/1e6,ray_datenum)
    odens[ci] = idens[0][0]
    hdens[ci] = idens[0][1]

crs_ll = convert2(crs_l, tvec_datetime,'GEO','car',['m','m','m'],'GEO','sph',['m','m','m'])
lats = [crl[1] for crl in crs_ll]

plt.plot(lats, edens)
plt.grid()
plt.xlabel('mlat')
plt.ylim([6, 12])
plt.show()
plt.close()
"""
# single testing
ray_start = [5824677.7777794544,       -6758550.8643073505,        4556168.2153918259]
#ray_start = convert2([ray_start],[ray_datenum],'SM','car',['m','m','m'],'GEO','car',['m','m','m'])

edens,idens = get_edens(np.array(ray_start)/1e6,ray_datenum)
print(10**(edens))

# confirmed that results are identical with the exception of zbrat....
# trust the fortran tracer more....