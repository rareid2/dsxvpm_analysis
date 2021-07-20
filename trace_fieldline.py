# austins thesis code
import numpy as np
import spacepy.irbempy as irbem
import spacepy.time as spt
import spacepy.coordinates as spc
from scipy.integrate import ode

# Bmodel index number (from irbempy)
#             - 0 = IGRF
#             - 1 = Eccentric tilted dipole
#             - 2 = Jensen&Cain 1960
#             - 3 = GSFC 12/66 updated to 1970
#             - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
#             - 5 = Centred dipole

def B_dir_3D(t, x, dtime, bmodel, extfield, direction):
    # p0 must be in GEO car in RE
    # dtime must be a datetime object

    pos = spc.Coords([x[0],x[1], x[2]],'GEO','car')
    tv = spt.Ticktock(dtime)
    B = irbem.get_Bfield(tv, pos, extMag=extfield,options=[1,0,0,0,bmodel], omnivals=None) 
    
    Bmags = direction*B['Bvec']/B['Blocal']
    
    return [Bmags[0][0], Bmags[0][1], Bmags[0][2]]

def trace_fieldline_ODE_3D(p0, stopalt, dtime, bmodel, extfield, direction):
    R_E = 6372.1 # km

    x = []
    y = []
    z = []
    dt = 0.05 # step size 
    r = ode(B_dir_3D)
    r.set_integrator('vode')
    
    r.set_initial_value(p0,0)
    r.set_f_params(dtime, bmodel, extfield, direction)
    counts = 0
    while r.successful():
        r.integrate(r.t + dt)
        x.append(r.y[0])
        y.append(r.y[1])
        z.append(r.y[2])

        counts+=1

        # stop conditions
        if np.linalg.norm(r.y) < 1:
            print('hit the Earth')
            break
            
        if counts > 1000:
            print('max count')
            break
        
        if np.linalg.norm(r.y) < ((stopalt+R_E)/R_E):
            print('hit stopalt')
            break

    #return np.array(x), np.array(y), np.array(z)
    return x[-1], y[-1], z[-1]

def get_bfield_irbem(x,dtime,bmodel,extfield,direction):

    pos = spc.Coords([x[0],x[1], x[2]],'GEO','car')
    tv = spt.Ticktock(dtime)
    B = irbem.get_Bfield(tv, pos, extMag=extfield,options=[1,0,0,0,bmodel], omnivals=None) 
    
    Bmag = B['Blocal']
    
    return Bmag 