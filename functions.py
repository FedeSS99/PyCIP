# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 16:39:57 2022
Updated on Wed Feb 22 12:09:00 2023

@author: Elvia Patrcia Barrón Cano
"""
import numpy as np
from math import pi,sin,cos,sqrt

"""
   SMOOTH function implementation based on IDL
   https://www.l3harrisgeospatial.com/docs/SMOOTH.html
"""

def smooth(a,w):
    n = len(a)
    min = (w-1)//2
    max = n - (w+1)//2
    r = np.zeros(n)
    
    for i in range(n):
        if (i < min or i > max):
            if i == 0:
                r[i] = (2*a[0]+a[1])/3
            elif i == n-1:
                r[i] = (a[n-2]+2*a[n-1])/3
            else:
                r[i] = a[i]
        else:
            for j in range(w):
                r[i] += a[i+j-w//2]
            r[i] = r[i]/w
    
    return r

"""
    I M P L E M E N T A T I O N   OF   THE   S O L A R   W I N D   M O D E L
    DOCTORAL THESIS
    Análisis de observaciones del MEXART: bases para estudios de Centelleo Interplanetario
    Julio César Mejía Ambriz 2012
    Appendix B.3.4 Explicit Expression and Computational Methodology  
    ecn (B.16), p.176 
    
    https://ru.dgb.unam.mx/handle/DGB_UNAM/TES01000692911
"""

# def myfun(x0,f,pf,ar,width,elong):       # Enable in case of 2 variables
def myfun(x0,f,pf,width,elong):          # Enable in case of 3 variables
# def myfun(x0,f,pf,elong):                # Enable in case of 4 variables

    # Parameters to adjust
    """ Enable in case of 2 variables """
    ve = x0[0]                           # Solar Wind Speed    
    alpha = x0[1]                        # Alpha = Spectral Index
    """ Enable in case of 3 variables """
    ar = x0[2]                           # AR = Axial Ratio
    """ Enable in case of 4 variables """
    # width = x0[3]                        # Theta0 = Source Angular Width
                    
    
    au = 1.49597e11                      # An Astronomical Unit in meters
    e = pi/180*elong                     # Elongation in Radians
    p1 = au*sin(e)                       # p = Minimum distance between the Sun and the disturbed Solar Wind
    qy = (np.arange(500)+1)*1.0e-6       # Component Y of the three-dimensional Wave Number of the Irregularities [1e-6,5e-4]
    z = np.flip(-np.arange(41)*0.05*au)  # 41*0.05= 2.05 [-2.05au, 0]
    z = z + au*cos(e)                    # Z Axis = Line of sight

    l = 2.14823                          # Wavelength
    sizef = len(f)                       # Amount of data to fit
    Pzf = np.zeros([sizef,41])           # Temporary array to save the result of the integral in Qy
    PFm = np.zeros(sizef)                # Temporary array to save the result of the integral in Z
     
    for i in range(41):
        for j in range(sizef):
            # Integrating in Qy
            R = sqrt((z[i]**2.0+p1**2))  # Ecn (B.4) Heliocentric distance to the layer
            Vx = ve*p1/R			     # Ecn (B.3) Velocity in X
            qx = 2*pi*f[j]/Vx            # Ecn (B.2) Component X of the three-dimensional Wave Number of the Irregularities
            q2 = qx**2 + qy**2	         # Ecn q2=q^2 = qx^2 + qy^2 + (qz^2=0)
            q2i = qx**2 + qy**2/ar**2    # For Non-Isotropic Medium
            z0 = z + au*cos(e)           # Defining z0
            # Ecn. B16
            Pzf[j,i] = np.trapz(R**(-3)*
                                np.sin((q2*z0[i]*l)/(4.0*pi))**2 * 
                                np.exp(-q2 * ( z0[i]*width*(4.845e-6)/2.35 )**2 ) * 
                                (q2i)**(-alpha/2), qy)
     
    for j in range(sizef):
        # Integrating in Z
        PFm[j] =  np.trapz(Pzf[j],z)

    PFm = PFm/PFm[0]                     # Normalizing

    return PFm - pf                      # Returning the difference between the modeled and the real values

