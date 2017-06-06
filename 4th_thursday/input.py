import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rc
from pylab import *
from scipy import *
from scipy.optimize import leastsq
from scipy.optimize import fmin_slsqp
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.optimize import newton
import os.path, os, sys
plt.rcParams.update({'font.size': 16, 'figure.figsize': [8, 6], 'legend.fancybox': True,
                     'xtick.major.size':6, 'ytick.major.size':6, 'figure.dpi':600,
                     'xtick.direction':'out','ytick.direction':'out'})
workdir=%pwd
print(workdir)







# Coordinates
x1 =[]
y1 =[]
x2 =[]
y2 =[]

# Starting coordinates
x1.append(4.) 
y1.append(5.)
x2.append(0.) 
y2.append(5.)

# Coordinates (isolated particles)
xS1 =[]
yS1 =[]
xS2 =[]
yS2 =[]

# Starting coordinates (isolated particles) = same as before
xS1.append(x1[0]) 
yS1.append(y1[0])
xS2.append(x2[0]) 
yS2.append(y2[0])


# Velocities
vx1 = []
vy1 = []
vx2 = []
vy2 = []

# Counter
i=0
# Viscosity
eta = 1.0
# Time step
dt = 0.2
# Time counter
iter = 20000
# Time = (Time counter)*(Time step)

# Particle 1
F1x = 0.
F1y = -1.
a1 = 1 
p1a = 1./(6.*pi*eta*a1)

# Particle 2
F2x = 0.
F2y = -1.0
a2 = 1
p2a = 1./(6.*pi*eta*a2)

# Off-diagonal element -> point-force (stokeslet)
pD = 1./(8.*pi*eta)

# Final coordinates of the isolated particles
xS1.append(xS1[0]+iter*dt*F1x*p1a) 
yS1.append(yS1[0]+iter*dt*F1y*p1a)
xS2.append(xS2[0]+iter*dt*F2x*p2a) 
yS2.append(yS2[0]+iter*dt*F2y*p2a)

for i in range(0,iter):
   
    # Relative distance between particle 1 and 2 in x and y-directions
    dx12 = x2[i]-x1[i]
    dy12 = y2[i]-y1[i]
    r12 = sqrt(dx12*dx12+dy12*dy12)
    
    # Matrix element
    # xx,xy, and yy 
    ExxD12 = pD*(1.+(dx12*dx12/(r12*r12)))/r12
    ExyD12 = pD*(dx12*dy12/(r12*r12))/r12
    EyyD12 = pD*(1.+(dy12*dy12/(r12*r12)))/r12

    Exx1 = p1a
    Exy1 = 0.
    Eyy1 = p1a
    
    Exx2 = p2a
    Exy2 = 0.
    Eyy2 = p2a
    
    vx1a = (Exx1*F1x+Exy1*F1x)+(ExxD12*F2x+ExyD12*F2y)
    vy1a = (Eyy1*F1y+Exy1*F1y)+(EyyD12*F2y+ExyD12*F2x)

    vx2a = (Exx2*F2x+Exy2*F2x)+(ExxD12*F1x+ExyD12*F1y)
    vy2a = (Eyy2*F2y+Exy2*F2y)+(EyyD12*F1y+ExyD12*F1x)
    

    vx1.append(vx1a)
    vy1.append(vy1a)
    vx2.append(vx2a)
    vy2.append(vy2a)
    
    x1.append(x1[i]+vx1[i]*dt)
    y1.append(y1[i]+vy1[i]*dt)
    x2.append(x2[i]+vx2[i]*dt)
    y2.append(y2[i]+vy2[i]*dt)
    

i = iter-1
vx1.append(vx1[i])
vy1.append(vy1[i])
vx2.append(vx2[i])
vy2.append(vy2[i])

x1 = asarray(x1)
y1 = asarray(y1)
vx1 = asarray(vx1)
vy1 = asarray(vy1)

x2 = asarray(x2)
y2 = asarray(y2)
vx2 = asarray(vx2)
vy2 = asarray(vy2)

# Scalar to rescale vectors (ss) and points each (ee)
ee = 100
ss = 10.

plt.quiver(x1[::ee],y1[::ee],vx1[::ee],vy1[::ee],scale=ss)
plt.quiver(x2[::ee],y2[::ee],vx2[::ee],vy2[::ee],scale=ss,color='r')
plt.plot(xS1,yS1,ls=':',color='k')
plt.plot(xS2,yS2,ls=':',color='r')

plt.ylabel(r"$\Delta x$",fontsize='12')
plt.xlabel(r"$\Delta y$",fontsize='12')

plt.subplots_adjust(left=0.20)
plt.subplots_adjust(right=0.96)
plt.subplots_adjust(bottom=0.14)
plt.subplots_adjust(top=0.90)
plt.savefig(workdir+'/Flow_sphere_lab.pdf')
plt.show()
