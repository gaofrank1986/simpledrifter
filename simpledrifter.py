#!/usr/bin/env python

import numpy as np

### Constants

g = 9.8  # Gravitational acceleration [m/s^2]

###

class Particle(object):
    """
    A water particle class.
    """
    def __init__(self,x0=0,z0=0,t0=0):
        """
        Particle instance constructor. x0, z0, t0 are initial positions
        in horizontal and vertical space and time.
        """
        self.x = np.array([x0]) # x-position
        self.z = np.array([z0]) # z-position
        self.t = np.array([t0]) # time

    def advect(self,u,w,dt):
        """
        Advects the particle forward in time with velocity (u,w) and 
        time step dt.
        """
        self.x = np.append(self.x,self.x[-1]+u*dt)
        self.z = np.append(self.z,self.z[-1]+w*dt)
        self.t = np.append(self.t,self.t[-1]+dt)


class Drifter(object):
    """
    A finite size drifter class.
    """ 
    def __init__(self,x0=0,z0=0,t0=0,H=1,A=1,Cd=1,m=0.5,tilt0=0,\
                 enableTilting=False):
        """
        Drifter instance constructor. 
        """ 
        self.x = np.array([x0])       # x-position
        self.z = np.array([z0])       # z-position
        self.t = np.array([t0])       # time
        self.tilt = np.array([tilt0]) # tilt

        self.Cd = Cd # drag coefficient
        self.A  = A  # surface area
        self.m  = m  # mass
        self.H  = H  # drogue height
        self.enableTilting = enableTilting
 
        # Drifter velocity
        self.u = np.array([0])
        self.w = np.array([0])

        # Angular velocity
        self.angvel = np.array([0])

        # Force acting on the drifter
        self.F_x = np.array([0])
        self.F_z = np.array([0])

    def advect(self,u,w,dt,a,k,rho=1E3):
        """
        Advects the object forward in time with velocity (u,w) and 
        time step dt.
        """
        t0 = self.tilt[-1]

        omega = np.sqrt(g*k)

        # Flow velocity relative to the object
        urel,wrel = u-self.u[-1],w-self.w[-1]

        F_x = np.sign(urel)*0.5*rho*urel**2*self.Cd*self.A*np.cos(t0)
        F_z = np.sign(wrel)*0.5*rho*wrel**2*self.Cd*self.A*np.sin(t0)

        self.F_x = np.append(self.F_x, F_x*np.cos(t0)-F_z*np.sin(t0))
        self.F_z = np.append(self.F_z,-F_x*np.sin(t0)+F_z*np.cos(t0))

        if self.enableTilting:

            # Tangential force around the pivot point
            Ft = F_x*np.cos(t0)-F_z*np.sin(t0)

            # Angular acceleration and velocity
            angacc = Ft/(0.5*self.m*self.H)
            angvel = self.angvel[-1] + angacc*dt       

            # Compute tilt based on angular velocity and acceleration
            t1 = t0 + angvel*dt + 0.5*angacc*dt**2

            self.tilt = np.append(self.tilt,t1)
            self.angvel = np.append(self.angvel,angvel)

        # Drogue center acceleration
        ax,az = self.F_x[-1]/self.m,self.F_z[-1]/self.m

        # Velocity of the drogue
        self.u = np.append(self.u,self.u[-1]+ax*dt)
        self.w = np.append(self.w,self.w[-1]+az*dt)

        u0,w0 = self.u[-1],self.w[-1]

        # Drogue center position
        x1 = self.x[-1] + u0*dt + 0.5*ax*dt**2
        #z1 = self.z[-1] + w0*dt + 0.5*az*dt**2
       
        # Force drifter to follow the surface
        t = self.t[-1]+dt
        eta = a*np.cos(k*x1-omega*t)
        z1 = eta-0.5*self.H

        # Update position and time
        self.x = np.append(self.x,x1)
        self.z = np.append(self.z,z1)
        self.t = np.append(self.t,self.t[-1]+dt)
