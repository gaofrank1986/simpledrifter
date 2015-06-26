#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from simpledrifter import Particle,Drifter

### Constants

g   = 9.8  # Gravitational acceleration [m/s^2]
rho = 1E3  # Water density [kg/m^3]

### Wave parameters

a = 0.2 # Wave amplitude [m]
k = 1   # Wavenumer [rad/m]

omega = np.sqrt(g*k) # Angular frequency [rad/s]
cp = omega/k         # Phase speed [m/s]
T = 2*np.pi/omega    # Period [s]

### Simulation parameters

dt       = 1E-3 # Time step [s]
duration = 30   # Simulation run time [s]

x     = np.linspace(0,2*np.pi,1E5) 
times = np.arange(0,duration+dt,dt) 

outputInterval = 1E2   # [s]

### Drifter parameters

H  = 0.6  # Drogue height [m]
A  = 0.36 # Drogue surface area [m^2]
Cd = 2    # Drag coefficient
m  = 0.5  # Drifter mass [kg]

t0 = 0                             # Initial deploy time
x0 = 0.5*np.pi                     # Initial position in x
z0 = a*np.cos(k*x0-omega*t0)-0.5*H # Initial position in z
tilt0 = 0.0                        # Initial tilt [rad]

enableTilting = False

###

# Initialize water particles to track
particles = []
for z_release in np.arange(-2.75,-0.,0.25):
    particles.append(Particle(x0=0.5*np.pi,z0=z_release,t0=0))

# Initialize drifter
d = Drifter(x0=x0,z0=z0,t0=t0,H=H,A=A,Cd=Cd,m=m,tilt0=tilt0,\
            enableTilting=enableTilting)

Hinv = 1./H # Inverse drogue height [1/m]

outputTime = 0
writeImage = True

for t in times:

    if t >= outputTime:
        writeImage = True

    # Output current snapshot 
    if writeImage:

        print 'Outputting snapshot at t = '+'%05.2f' % t

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111,xlim=(0,2*np.pi),ylim=(-3,1))
        plt.plot(x,a*np.cos(k*x-omega*t),'k-')
        for p in particles:
            plt.plot(p.x,p.z,'r-',lw=0.2)
            plt.plot(p.x[-1],p.z[-1],'r.')
        plt.plot([p.x[-1] for p in particles],\
                 [p.z[-1] for p in particles],'b:')
        plt.plot(d.x,d.z,'k-',lw=0.4)
        plt.plot(d.x[-1],d.z[-1],'ko')
        plt.plot([d.x[-1]-0.5*H*np.sin(d.tilt[-1]),\
                  d.x[-1]+0.5*H*np.sin(d.tilt[-1])],\
                 [d.z[-1]-0.5*H*np.cos(d.tilt[-1]),\
                  d.z[-1]+0.5*H*np.cos(d.tilt[-1])],\
                 'k-',lw=2)
        plt.grid(True)
        plt.xlabel('x [m]',fontsize=16)
        plt.ylabel(r'$\eta$ [m]',fontsize=16)
        plt.title('t = '+'%05.2f' % t+' s',fontsize=16)
        plt.savefig('eta_'+'%05.2f' % t+'.png',dpi=100)
        plt.close(fig)

        writeImage = False
        outputTime += outputInterval

    # Advect water particles forward
    for p in particles:
        xp = p.x[-1]
        zp = p.z[-1]
        u = a*omega*np.exp(k*zp)*np.cos(k*xp-omega*t)
        w = a*omega*np.exp(k*zp)*np.sin(k*xp-omega*t)
        p.advect(u,w,dt)

    # Water velocity integrated over drogue height
    u = a*cp*Hinv*(np.exp(k*d.z[-1]+0.5*H)-np.exp(k*d.z[-1]-0.5*H))\
                 *np.cos(k*d.x[-1]-omega*t)
    w = a*cp*Hinv*(np.exp(k*d.z[-1]+0.5*H)-np.exp(k*d.z[-1]-0.5*H))\
                 *np.sin(k*d.x[-1]-omega*t)

    # Advect the drifter forward
    d.advect(u,w,dt,a,k,rho=rho)

# Flow velocity integrated over drifter height
uw = a*cp*Hinv*(np.exp(k*d.z+0.5*H)-np.exp(k*d.z-0.5*H))\
              *np.cos(k*d.x-omega*d.t)

# Flow acceleration
ac = omega*a*cp*Hinv*(np.exp(k*d.z+0.5*H)-np.exp(k*d.z-0.5*H))\
              *np.sin(k*d.x-omega*d.t)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
plt.plot(d.t,uw,'r-',label='Water')
plt.plot(d.t,d.u,'k-',label='Drifter')
plt.plot(d.t,uw-d.u,'g-',label='Water-drifter')
plt.grid(True)
plt.xlabel('Time [s]',fontsize=16)
plt.ylabel('u-velocity [m/s]',fontsize=16)
plt.ylim(-0.6,0.6)
plt.xlim(0,6)
plt.legend(loc='lower left',fancybox=True,shadow=True,ncol=3)
plt.title('Water and drifter u-velocity',fontsize=16)
plt.savefig('uvel2.png',dpi=100)
plt.close(fig)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
plt.plot(d.t,ac,'r-',label='Water')
plt.plot(d.t,d.F_x/d.m,'k-',label='Drifter')
plt.plot(d.t,ac-d.F_x/d.m,'g-',label='Water-drifter')
plt.grid(True)
plt.xlabel('Time [s]',fontsize=16)
plt.ylabel('Acceleration [m/s^2]',fontsize=16)
plt.ylim(-2,2)
plt.xlim(0,6)
plt.legend(loc='lower left',fancybox=True,shadow=True,ncol=3)
plt.title('Water and drifter acceleration',fontsize=16)
plt.savefig('accel2.png',dpi=100)
plt.close(fig)

