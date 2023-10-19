# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 11:31:10 2023
@author: Elvia Patricia Barrón Cano

auto-mexart-20230217.pro
@autor: Julio Cesar Mejía Ambriz
last  20230217

This program gets:
    
1. Parameters for the WIPSS FORMAT
2. ROW and SMOOTH SPECTRUM DATA
3. VELOCITY MODEL CHARTS
4. IPS VELOCITY MAP

For the NEW MEXART DATA (2020 onwards)
"""
import os,sys
import numpy as np
from datetime import datetime
from time import ctime,gmtime,localtime,strftime
from functions import smooth,myfun
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from waveletFunctions import wavelet
from math import pi,sin,cos,asin,acos,tan,sqrt,isnan
import scipy.optimize as opt

"""
  The current directory which is where the running program resides,
  is denoted by "." 
 
   At that level is the wipss subdirectory.
 
  And in wipss they are: data, list, spec-data-smooth, spectra-data & 
  spectra-model. Where the data to be analyzed is located and the reports
  and graphs will be stored.
"""

t0 = ctime()
print('\n     A  U  T  O    -    M  E  X  A  R  T\n')
print('BEGINING at %s' %(t0))
path =  './wipss/'
data_dir = path + 'logs/'
i = 0

if not os.path.exists(data_dir):
    print('Making Logs Dir ... %s' %(data_dir))
    os.makedirs(data_dir)
    i += 1

mylocaltime = strftime("%Y%m%d-%H%M%S",localtime())  
fn1 = data_dir+'log-'+mylocaltime+'.txt'
o0 = open(fn1,'w')              # Opening the Output log file for writing

o0.write('\n     A  U  T  O    -    M  E  X  A  R  T\n\n')
o0.write('BEGING at %s\n' %(t0))
if (i == 1):
    o0.write('Making Logs Dir: %s\n' %(data_dir))
        
# The Absolute Path for DATA ACCESS is defined
data_dir = path + 'data/'

"""
  A FOR is performed on the objects contained in the directory
 
  The OS.LISTDIR function gives us a complete list of files and folders
  that they live in that directory. Only the names, without their path.
  
  It is validated that only those that are FILES and that end in .dat
  are added to the list, with the OS.PATH.ISFILE functions and the 
  *.ENDSWITH('.dat') property.
"""

files = sys.argv
files.remove(sys.argv[0])
nf = len(files)

# files = [fn1 for fn1 in os.listdir(data_dir) if (os.path.isfile(data_dir+fn1) & fn1.endswith('.dat'))]
# #files = [fn1 for fn1 in os.listdir(data_dir) if (os.path.isfile(data_dir+fn1) & fn1.endswith('.txt'))]
# nf = len(files)         # number of files = number of sources

"""
 Defining Position Names of SUN and SOURCES
"""

# Sun
raSun = np.zeros(nf)          # solar right ascension (AR) at time of observation
la = np.zeros(nf)             # solar latitude at the time of observation

# Source
raSource = np.zeros(nf)       # right ascensions
delta = np.zeros(nf)          # declinations

# Sun + Source
elong = np.zeros(nf)          # elongations
arc = np.zeros(nf)            # to plot angles
he = np.zeros(nf)             # heliocentric latitudes at the time of observation
res = np.zeros(nf)            # difference between RA of SOURCE and RA of SUN 
S = np.empty(nf,dtype='U1')   # West (W) or East(E) side
p_angle = np.zeros(nf)        # Position Angle from solar North (same as LASCO)

Tres = np.zeros((nf,365))     #
Telong = np.zeros((nf,365))   # elongations along the year
The = np.zeros((nf,365))      # heliocentric latitudes along the year
# TraSun = np.zeros(365)      # Solar right ascencions along the year. No need to define it

Tl = np.zeros(365)            # All Solar Latitudes

# conversion factors
r = pi/180.0                  # degrees => radians
d = 1/r                       # radians => degrees
# Days of the month
days_of_m = [31,28,31,30,31,30,31,31,30,31,30,31]
Source = np.empty(10,dtype='U1') # Radio Source Name
rsnp = 0                    # Number of Radio Sources NOT Analyzed
rsnan = 0                   # Number of Radio Sources with NAN
rsnid = 0                   # Number of Radio Sources NOT Identified

""" 
   MAIN CYCLE
   
   In this cycle we calculate:
       All Angles and Distances, Time Series, 
       Wavelet, M-Index, Power Spectrums, 
       Velocities & Theoretical Spectrums
"""
print('THERE ARE %d files to process' %(nf))
o0.write(('THERE ARE %d files to process\n' %(nf)))
for m in range(nf):

    print('\nPROCESSING the file #%d/%d ... %s' %(m+1,nf,files[m]))
    print('(1)  Computing elengation & others')
    o0.write('\nPROCESSING the file #%d/%d ... %s\n' %(m+1,nf,files[m]))
    o0.write('(1)  Computing elengation & others\n')
 
    #=======COMPUTING ELENGATION=============        
    year = int(files[m][0:4])       # year
    month = int(files[m][4:6])      # month
    di  = int(files[m][6:8])        # day of the mont
    date = files[m][0:8]            # date
    hh  = int(files[m][9:11])       # hour of RA
    mm  = int(files[m][11:13])      # min of RA
    ss  = int(files[m][13:15])      # sec of RA
    sign = files[m][15]             # sign
    dd  = int(files[m][16:18])      # dec in degrees
    mi  = int(files[m][18:20])      # dec in min
    sec = int(files[m][20:22])      # dec in sec
    if (sign == '-'):
        dd = -dd
        mi = -mi
        sec = -sec

    doy = 0   # Days accumulated until the beginning of the month. Where for jan 01 has a doy = 0
    for i in range(month-1) :
        doy= days_of_m[i]+doy
        # Bisiesto
        if (not(year % 4) and (i == 1) and ((year%100) or (not(year%400)))):
            doy += 1
            
    doy += di
    # Spencer's solar latitude in RADIANS Ecn. (2.7)
    #####################################################################################
    # DIFERENCIA CLAVE en elongation.pro por 1, y varia enlong, la, he, arc y p_angle
    #####################################################################################

    g =((2*pi)/365)*(doy-1)                  # Gamma diary angle in RADIANS
    # latitude  in RADIANS
    lat = (0.006918 - 0.399912*cos(g)+0.070257*sin(g)-0.006758*cos(2*g)+0.000907*sin(2*g)-0.002697*cos(3*g)+0.00148*sin(3*g))
    # latitude in DEGREE
    la[m] = lat*d
    # Right Ascension in DEGREES (24hrs=360° => 1hr= 15°; 60m=15° => 1m=0.25°; 60s=0.25° => 1s=0.00416667)
    raSource[m] = (hh*15)+(mm*0.25)+(ss*0.00416667)
    # right ascension in RADIANS
    alfa = raSource[m]*r
    # declination in DEGREES
    delta[m]=dd+mi/60+sec/3600
    # declination in RADIANS
    deltar = delta[m]*r
    # Solar Right Ascension in DEGREES. The vernal equinox started in march 20                        
    if ((doy >= 80) and (doy <= 365)) :
        raSun[m] = (doy-79)*0.9856
    else :
        raSun[m] = (doy+286)*0.9856
    # Solar Right Ascension in RADIANS
    theta = raSun[m]*pi/180
    # Ecn. (2.3) pag. 49                  
    e = cos(lat)*cos(deltar)*cos(theta-alfa)+sin(lat)*sin(deltar)
    # Elongation in RADIANS
    elongr = acos(e)
    # Elongation in DEGREES
    elong[m]=elongr*d

    if ((elong[m] > 18) and (elong[m] < 90)) :       
        # Heliocentric Angle in DEGREES, Ecn. (2.6) 
        he[m] = asin(sin(deltar-lat)/tan(elongr))*d
    
        res[m]=raSun[m]-raSource[m]
        if (abs(res[m]) <= 180): 
            if res[m] >0 : S[m]='W' 
            else: S[m]='E'
        else:
            if res[m] >0 :  S[m]='E'
            else: S[m]='W'
            res[m] = 360 - abs(res[m])
        
        #=================== ALL SOLAR LATITUDES =================================
        # Solar Latitudes Ecn. (2.7) pag. 49 
        Tg = np.arange(365)*(2*pi/365)
        # in RADIANDS
        Tlr = 0.006918 - 0.399912*np.cos(Tg)+0.070257*np.sin(Tg)-0.006758*np.cos(2*Tg)+0.000907*np.sin(2*Tg)-0.002697*np.cos(3*Tg)+0.00148*np.sin(3*Tg)
        # in DEGREES
        Tl = Tlr*d
        # Solar Right Ascensions along the year
        TraSun = np.concatenate(((np.arange(79)+287)*0.985647,(np.arange(286)+1)*0.985647))
        for i in range(365):
            Tres[m,i] = abs(TraSun[i] - raSource[m])
            if Tres[m,i] > 180: Tres[m,i] = 360 - Tres[m,i]
        # All Elongations Ecn. (2.3) pag. 49  
        Telong[m] = d*np.arccos(np.cos(Tlr)*np.cos(deltar)*np.cos(Tres[m]*r) + np.sin(Tlr)*np.sin(deltar))
        # All Heliocentric Angles
        for i in range(365):
            if Telong[m,i] <= 90:
                The[m,i] = np.arcsin(sin(deltar-Tlr[i])/tan(Telong[m,i]*r))*d
        # Calculate PA of P-point (deg) counterclockwise from North
        arcx = he[m]/(np.absolute(The[m]).max())*90
        if elong[m] > 90 :
            arc[m] = 999
        elif he[m] >= 0 : 
            if S[m] == 'W' : 
                arc[m] = arcx
            else :
                arc[m] = 180 - arcx
        else :
            if S[m] == 'W' : 
                arc[m] = 360 + arcx
            else :
                arc[m] = 180 - arcx
           
        p_angle[m] =  arc[m]
        if ((arc[m] < 0) or (arc[m] > 90)) :
            p_angle[m] -= 90 

        #################################################################################
        # OBTAINING TIME SERIES WITH 10s DETRENDING (aa vs times)
        #################################################################################
    
        print('(2)  Obtaining Time Series')
        o0.write('(2)  Obtaining Time Series\n')
        i1 = np.loadtxt(data_dir+files[m], dtype={'names':('times','voltage'),
                              'formats':('U15','float')},usecols= (1,2))
        n = i1.size
        times = []
        for i in range(n):
            tt = datetime.strptime(i1[i]['times'],"%H:%M:%S.%f").time() # Converting the string to numbers in a tuple
            times.append(tt.hour*3600.0 + tt.minute*60.0 + (tt.second + tt.microsecond/1000000)) # Storing time in seconds
    
        voltage_wn = np.asarray([ x for x in i1['voltage'] if isnan(x) == False])
        for i in range(n):    
            if (isnan(i1[i]['voltage'])):
                i1[i]['voltage'] = voltage_wn.mean()
        
        nunan = n - voltage_wn.size         # Number of Nans
        print('     ==> There are(is) %d NANs' %(nunan))
        o0.write('     ==> There are(is) %d NANs\n' %(nunan))
        if (nunan != 0):
            rsnan += 1  

        source = files[m].removesuffix('.dat')
        # source = files[m].removesuffix('.txt')
        source = source.split(sep="_")[1]   # source radio name from filename
        n_time = len(times)                 # length of times list
        MidObsUT = strftime("%T ",gmtime(times[n_time//2]))  # 'hh:mm:ss' transit peak
        times = np.asarray(times)           # converting from list to array
        # DT = the sampling time = 50 samples per second
        dt = 0.02                           # (2.38419e-07)*3600*24
        aa = i1[:]['voltage'] - smooth(i1[:]['voltage'],500) # 10 sec filter
        # plt.plot(voltage)
        # plt.title('Voltage')
        # plt.show()
        # plt.plot(aa)
        # plt.title('La diferencia = aa')
        # plt.show()
        """
          OBTAINING WAVELET
            wave = WAVELET(Y,DT)
          INPUTS:
            Y = the time series of length N.
            DT = amount of time between each Y value, i.e. the sampling time.
        """
        print('(3)  Obtaining Wavelet')
        o0.write('(3)  Obtaining Wavelet\n')
        pad = 1         # pad => PAD
        s0 =  2*dt      # shortest scale   
        dj = 0.0625     # spacing between discrete scales. Default is 0.125  
        j1 = 14/dj      # of scales minus one# Default is j = alog2(n*(dt/s0))/dj 
        mother = 'MORLET'
        """
          Note: for accurate reconstruction and variance computation, set: 
            s0 = dt    for Morlet
            s0 = dt/4  for Paul
            (Most commonly, s0=2*dt
          wave = WAVELET(aa,dt,PERIOD=period,SCALE=scale,S0=s0,COI=coi,DJ=dj,J=j1,MOTHER=mother,/RECON,/PAD,signif=signif)
        """
        wave, period, scale, coi = wavelet(aa,dt,pad,dj,s0,j1,mother)   # FALTA /RECON,signif=signif
        power = (np.abs(wave))**2                                       # compute Wavelet POWER SPECTRUM
        # plt.plot(period,wave)
        # plt.plot(power)
        # plt.title('Wavelet Power Spectrum')
        # plt.show()
        #=========================================================================

        global_ws = (np.sum(power, axis=1) / n)  # global wavelet spectrum (GWS), time-average over all times
 
        # We define BOX, in time, one minute before and one minute after the EMISSION PEAK.
        dt1    = times.max()-times.min()
        dt1med = dt1/2
        ax1    = times.min()+dt1med-60
        ax2    = times.min()+dt1med+60
        #===============================================================
        # M-INDEX
        #===============================================================
        print('(4)  Calculating the M-Index')
        o0.write('(4)  Calculating the M-Index\n')
        ww = np.argwhere((times >= ax1) & (times <= ax2))
        ww2 = times[ww]
        pp= period[np.where((period >= 1) & (period <= 3.3))]
        a = power[73:102,ww[0,0]:ww[-1,0]+1]        # 1 Hz a 0.3 Hz
        aoff = power[20:58,ww[0,0]:ww[-1,0]+1]      # 10.17 Hz a 2 Hz
        # aall = power[73:127,ww[0,0]:ww[-1,0]+1]     # 1.02417 Hz a 0.103088 Hz. No used
        mvalue = np.sqrt(np.mean(a)-np.mean(aoff))  # CALCULO DEL INDICE PARA WIPSS
        # plt.plot(a)
        # plt.title('a')
        # plt.show()
        # plt.plot(aoff)
        # plt.title('aoff')
        # plt.show()
        # plt.plot(aall)
        # plt.title('aall')
        # plt.show()
    
        #==== Generating files with information from POWER SPECTRUM DATA & SMOOTH ====
        
        print('(5)  Generating the Power Spectrum Data file')
        o0.write('(5)  Generating the Power Spectrum Data file\n')
        fn1 = path+'spectra-data/'+source+'-'+date+'-spectrum.dat'
        o1 = open(fn1, 'w')              # Opening the Output file for writing
        a2 = power[20:127,ww[:,0]]
        powr = np.zeros(107)
        frq = np.zeros(107)
                        
        for i in range(107) :
            powr[i] = np.mean(a2[i])
            frq[i] = 1/period[20+i]
            # Saving: Frecuency & Power 
            o1.write('%15.10f  %25.5f\n' % (frq[i],powr[i]) )
        o1.close()                       # Closing Output file 

        print('(6)  Generating the Power Spectrum Data Smooth file')
        o0.write('(6)  Generating the Power Spectrum Data Smooth file\n')
        powr2 = smooth(np.flip(powr),10) # Power Smoothing
        frq2 = np.flip(frq)
        i = frq2[frq2 < 0.302504].size
        X = frq2[i:i+13]
        Y = powr2[i:i+13]
        AJUSTE_BF = np.polyfit(X, Y, 1)  # Linear Fit
        p2 = np.zeros(24) 
        for i in np.arange(23,-1,-1) :
            p2[i] = AJUSTE_BF[0]*frq2[i] + AJUSTE_BF[1]

        powr2[0:24] = p2
        fn1 = path+'spec-data-smooth/'+source+'-'+date+'-spectrum-smooth.dat'            
        o1 = open(fn1, 'w')              # Opening the Output file for writing
        for i in range(107) :
            # Saving: Frecuency & Power 
            o1.write('%15.10f  %25.5f\n' % (frq2[i],powr2[i]) )
                
        o1.close()                       # Closing Output file 
        # plt.plot(frq2,powr2)
        # plt.show()
        # FINISH part of WAVELETS
        
        #=============================================================================
        #                         CALCULATING SPEED
        #=============================================================================
        
        print('(7)  Calculating SPEED\n')
        o0.write('(7)  Calculating SPEED\n')
        l = 2.14823                             # WAVELENGTH for MEXART
        instrument ='MEXART'                    # Instrument Name

        """
            These ones are the initial values (ar=1.5 comment by coles)
            velocity in km/sec, alpha, ar, angular size & source elongation (this only one no change)
        """

        x0 = [285000,3.5,1]                      # Initial Conditions: Velocity, Alpha & Axial Ratio
        # x0 = [200000, 1.5]
        # x0 = [285000, 2.6, 1.0, 0.25]
        # x0 = [300000, 3.3, 1.0]
        # ar = 1.0                               # AR = Axial Ratio
        width = 0.2                              # Theta0 = Source Angular Width
        ex = elong[m]                            # Elongation
     
        limits = ([100000,3.3,0.7],[2000000,3.8,1.3])    # Limits
        # limits = ([200000,1.5,-5,0.001],[1500000,6.8,10,0.5])
        #limits = ([200000,1.5,-5,0.05],[1500000,6.8,5,0.4])       
        intervalo = np.argwhere((frq2 > 0.1) & (frq2 < 1.2))
        intervalo = intervalo.reshape(intervalo.size)
        f = frq2[intervalo]                      # limiting the Frequency
        pf = powr2[intervalo]                    # and the power in that range
        pf = smooth(pf,5)                        # Power Smoothing
        pf = pf/pf[0]                            # Normalizing
        pf2 = powr2                              # A COPY WITHOUT limits
        pf2 = smooth(pf2,5)                      # Power Smoothing
        pf2 = pf/pf[0]                           # Normalizing
        """
          Theoretical Calculation of the SOLAR WIND (ECN 5.22, P.115) 
        
          FIT using LEAST SQUARES with the function MYFUN
          
          Method ‘lm’ : Levenberg-Marquardt algorithm. NOTE: This one doesn't support limits.
          Method ‘trf’ : Trust Region Reflective algorithm. NOTE: This one support limits.
        """
        # Varying only: Velocity & Alpha
        # pfit = opt.least_squares(myfun,x0,method='lm',verbose=2,args=(f,pf,ar,width,ex))

        # pfit = opt.least_squares(myfun,x0,bounds=limits,verbose=2,args=(f,pf,ar,width,ex))

        # Varying only: Velocity, Alpha & AR
        pfit = opt.least_squares(myfun,x0,bounds=limits,verbose=2,args=(f,pf,width,ex))
        # Varying:  Velocity, Alpha, AR & Width
        # pfit = opt.least_squares(myfun,x0,bounds=limits,verbose=2,args=(f,pf,ex))
        
        o0.write('\n' + pfit.message + '\n')
        o0.write('Function evaluations %d, final cost= %e & first-order optimality= %e.\n' %(pfit.nfev,pfit.cost,pfit.optimality))
        
        #======================= RECONSTRUCTING the THEORETICAL SPECTRUM ===================
        print('\n(8)  Reconstructing the Theoretical Spectrum')
        o0.write('\n(8)  Reconstructing the Theoretical Spectrum\n')
        
        ve = pfit.x[0]                          # Solar Wind Speed
        alpha = pfit.x[1]                       # Alpha = Indice Espectral
        """ Enable in case of 3 variables """
        ar = pfit.x[2]                        # AR = Axial Ratio
        """ Enable in case of 4 variables """
        # width = pfit.x[3]                     # Theta0 = Source Angular Width

        e = ex*pi/180                           # Elongation in Radians
        l = 2.1482                              # Wavelength
        au = 1.49597e11                         # An Astronomical Unit in meters
        p = au*sin(e)                           # p = Minimum Distance between the SUN and the DISTURBED SOLAR WIND
        z = np.flip(-np.arange(41)*0.05*au)     # 41*0.05= 2.05 [-2.05au, 0]
        z = z + au*cos(e)                       # Z-axis = Line of Sight
        qy = (np.arange(500)+1)*1.0e-6          # Component Y of the three-dimensional Wave Number of the Irregularities [1e-6,5e-4]
 
        size2 = len(f)                          # Amount of data to fit
        Pzf = np.zeros([size2,41])              # Temporary array to save the result of the integral in Qy
        PFm = np.zeros(size2)                   # Temporary array to save the result of the integral in Z
     
        for i in range(41):
            for j in range(size2):
                # integrating in Qy
                R = sqrt((z[i]**2.0+p**2))      # Ecn (B.4) Heliocentric Distance to the layer
                Vx = ve*p/R	                    # Ecn (B.3) Velocity in X
                qx = 2*pi*f[j]/Vx               # Ecn (B.2) Component X of the three-dimensional Wave Number of the Irregularities
                q2 = qx**2 + qy**2	            # Ecn q2=q^2 = qx^2 + qy^2 + (qz^2=0)
                q2i = qx**2 + qy**2/ar**2       # for Non-Isotropic Medium
                z0 = z + au*cos(e)              # Defining z0

                Pzf[j,i] = np.trapz(R**(-3)*
                    np.sin((q2*z0[i]*l)/(4.0*pi))**2 * 
                    np.exp(-q2 * ( z0[i]*width*(4.845e-6)/2.35 )**2 ) * 
                    (q2i)**(-alpha/2), qy)
     
        for j in range(size2):
            # Integrating in Z
            PFm[j] =  np.trapz(Pzf[j],z)

        pfited  = PFm/PFm[0]                    # Normalizing
        
        #============================ CALCULATING CHI SQUARE ==========================
        print('(9)  Calculating Chi Square')
        o0.write('(9)  Calculating Chi Square\n')
        chisy = np.zeros(size2)
        for i in range(size2):
            chisy[i] = (pf[i]-pfited[i])**2 / pfited[i]
        
        chisq = chisy.sum()
        
        #==================== Plotting the REAL & THEORETICAL SPECTRUM ===================
        print('(10) Plotting the REAL & THEORETICAL SPECTRUM')
        o0.write('(10) Plotting the REAL & THEORETICAL SPECTRUM\n')
        e = elong[m]*pi/180                     # Elongation in Radianes
        p1 = sin(e)                             # Distance
        
        fn1 = './radio sources.txt'
        fn2 = path+'spectra-model/'+date+'-'+source+'.png'  # Output files (PNG,EPS)
        i1 = np.loadtxt(fn1, dtype={'names':('name','ra_dec'),'formats':('U10','U13')},usecols= (0,1),skiprows=1)
        i = np.argwhere(i1['ra_dec']== source)
        if ( len(i) == 0):
            Source[0] = 'J' 
            Source[1:5] = list(source)[0:4]
            Source[5:10] = list(source)[6:11]
            Source=''.join(Source)
            rsnid += 1
            print('     * Radiofuente NO localizada en la Base de Datos *')
            o0.write('     * Radiofuente NO localizada en la Base de Datos *\n')
        else:
            Source = i1['name'][i][0][0].center(10)
        fig = plt.figure(dpi=150)
        ax = fig.add_subplot(111)
        # plt.plot(fi,pf,'b',label='obs')
        plt.plot(f,pf2,'k-',label='obs')
        plt.plot(f,pfited,'k--',label='fit')
        ax.set_title('MEXART@139.65 MHz'+' '+Source+' '+date[0:4]+'/'+date[4:6]+'/'+date[6:8])
        plt.ylabel('Normalized power')
        plt.xlabel('Frequency [Hz]')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(pf.min(), 6)
        plt.xlim(0.1, 3)
        plt.legend(loc='upper left')
        # plt.text(0.5, 4.5, 'AR = %3.1f' %(ar))
        plt.text(0.5, 4.5, 'AR = %3.1f' %(pfit.x[2]))
        plt.text(1.0, 4.5, r'$\chi^2 = %f$' %(chisq))
        plt.text(0.5, 3.5, r'$\alpha = %5.3f$' %(pfit.x[1]))
        plt.text(1.0, 2.5, 'V = %8.3f km/s' %(pfit.x[0]*1e-3))
        plt.text(0.5, 2.5, r'$\theta = %4.2f$' %(width))
        # plt.text(0.5, 2.5, r'$\theta = %4.2f$' %(pfit.x[3]))
        plt.text(0.7, 1.5, r'$\epsilon = %5.3f°$' ' (%5.3f AU)' %(elong[m],p1))
        ax.set_rasterized(True)
        plt.savefig(fn2)                        # saving the graph as PNG with 150dpi(ebook)

        #============ Generating file with information for WIPSS FORMAT =========
        
        print('(11) Adding an Entry to the WIPSS FILE')
        o0.write('(11) Adding an Entry to the WIPSS FILE\n')
        fn2 =  path+'list/WIPSS-FORMAT-FILE-'+mylocaltime+'.dat'
        if (os.path.isfile(fn2)):    
            o1 = open(fn2, 'a')   # Opening to add DATA to the Output file
        else:
            o1 = open(fn2, 'w')   # Creating the file and opening it for writing
            o1.write('Date      MidObsUT  Dur. Site  Freq  BW     Source    Size  RA-J2000  Dec-J2000  Limb  Dist.   Lat.   PA     Elong   Vel.   V-err  g-value      g-err  Method     Vel.  V-err g-value g-err Method\n')
      
        # date_wipss = date     # No need to define it
        MidObsUT   = MidObsUT
        Dur        = '2.0'
        Site       = 'MXRT'
        Freq       = '140'
        BW         = 12.5
        Size       = width      # 0.2
        RA         = source[0:2]+' '+source[2:4]+' '+source[4:6]
        DC         = source[6:9]+' '+source[9:11]+' '+source[11:13]
        Limb       = S[m]
        Dist       = 215*sin(elong[m]*r)
        Lat        = he[m]
        PA         = p_angle[m]
        Elonga     = elong[m]
        Vel        = pfit.x[0]*1e-3
        Verr       = -999
        # mvalue     = mvalue   # No need to define it
        gerr       = -999
        Method     = 'SS'
        
        # Saving
        o1.write('%s  %s %s  %s  %s  %s  %s  %s   %s  %s   %s   %6.2f %6.2f  %6.2f %6.2f  %6.1f  %d %13.5e  %d   %s\n'
                 % (date,MidObsUT,Dur,Site,Freq,BW,Source,Size,RA,DC,
                    Limb,Dist,Lat,PA,Elonga,Vel,Verr,mvalue,gerr,Method))

        o1.close()  # Closing Output file
    
    else:
        print('\n  * The ELONGATION parameter is OUT OF RANGE (< 18° or > 90°)*')
        print('          Elongation = %d' %(elong[m]))
        print('      The file %s was not processed' %(files[m]))
        o0.write('\n  * The ELONGATION parameter is OUT OF RANGE (< 18° or > 90°)*\n')
        o0.write('          Elongation = %d\n' %(elong[m]))
        o0.write('      The file %s was not processed\n' %(files[m]))
        rsnp += 1
        
# END of MAIN CYCLE 

"""
  Plotting RADIO SOURCES MAP
"""

fn1 = path+'list/WIPSS-FORMAT-FILE.dat'   # Input Filename
if (os.path.isfile(fn1)): 
    print('\nPLOTTING RADIO SOURCES MAP ... ')
    o0.write('\nPLOTTING RADIO SOURCES MAP ... \n')
    
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }

    fn2 = path +'spectra-model/'+'maps-'+mylocaltime+'.png' # Output Filename

    i1 = np.loadtxt(fn1, dtype={'names':('Elonga','Vel'),
                    'formats':('float','float')},
                    usecols= (18,19), skiprows=1)

    if (i1.ndim == 0):                        # If it processes ONLY ONE FILE it RESIZES
        i1 = i1.reshape(1)
    m = i1.size                               # Number of entries in the file WIPSS 
    
    theta0 = np.arange(360)*r
    theta = np.arange(3600)*r

    r0 = np.zeros(360)
    # Radial Distances in Degrees
    r1 = np.zeros(3600)+20      # 20 deg 
    r2 = r1+20      # 40 deg
    r3 = r2+20      # 60 deg
    r4 = r3+20      # 80 deg
    r5 = r4+10      # 90 deg

    # Velocity & circle size 
    v = [(200,4),(300,5),(400,6),(500,7),(600,8),
         (700,9),(800,10),(900,11),(1000,12),(1100,13)]

    # Adding the legend elements
    legend_elements = []

    for i in range(len(v)-1):
        legend_elements.append(Line2D([0],[0], marker='o',color='w',label=v[i][0],markerfacecolor='b',markersize=v[i][1]))

    legend_elements.append(Line2D([0],[0], marker='o',color='w',label='>%d' %(v[9][0]),markerfacecolor='b',markersize=v[9][1]))

    fig = plt.figure(figsize=(9,5),dpi=150)     # Figure size
    ax = fig.add_subplot(111)                   # ONE Plot
    ax.legend(handles=legend_elements, loc='right') # Plotting Legend

    plt.axes(projection = 'polar',frameon=False)    # NEXT PLOT will be POLAR COORDINATES 
    ax.set_title('     Apparent position of sources in the sky',pad=20,fontdict=font)   # TITLE

    plt.polar(theta0,r0,'ro')                    # Circle with radius r0          
    plt.polar(theta,r1,'k',linewidth=1.5)       # Circle with radius r1  
    plt.polar(theta,r2,'k',linewidth=2)         # Circle with radius r2
    plt.polar(theta,r3,'k',linewidth=2)         # Circle with radius r3
    plt.polar(theta,r4,'k',linewidth=2)         # Circle with radius r4
    plt.polar(theta,r5,'k',linewidth=2)         # Circle with radius r5

    plt.figtext(0.5,0.59,'20°')                 # Text 20°
    plt.figtext(0.5,0.67,'40°')                 # Text 40°
    plt.figtext(0.5,0.75,'60°')                 # Text 60°
    plt.figtext(0.5,0.83,'80°')                 # Text 80°
    plt.figtext(0.5,0.87,'90°')                 # Text 90°
    plt.figtext(0.48,0.39,'0.34 AU')            # Text 0.34 AU
    plt.figtext(0.48,0.31,'0.64 AU')            # Text 0.64 AU
    plt.figtext(0.48,0.23,'0.86 AU')            # Text 0.86 AU
    plt.figtext(0.48,0.10,'1 AU')               # Text 1 AU

    # DATE
    plt.figtext(0.65,0.78,'%s-%s-%s' %(date[0:4],date[4:6],date[6:8]))
    # Legend Title
    plt.figtext(0.76,0.74,'IPS speed [km/s]')
    # According to the velocity will be the size of the circle
    for i in range(nf):
        size2 = v[0][1]
        for j in range(len(v)-1):
            if ((i1[i+m-nf]['Vel'] >= v[j][0]) and (i1[i+m-nf]['Vel'] < v[j+1][0])) : 
                size2 = v[j][1]
                break
            else:
                if i1[i+m-nf]['Vel'] >  v[9][0] :
                    size2= v[9][1]
                    break                   
            
        plt.polar(arc[i]*r,i1[i+m-nf]['Elonga'],'b',marker='o',markersize=size2)

    plt.xticks([])                          # removing xticks
    plt.yticks([])                          # removing yticks
    ax.set_axis_off()                       # removing the box 
    ax.set_rasterized(True)
    plt.savefig(fn2)                        # saving the graph as PNG with 150dpi(ebook)

else:
    print('\nNOTE: The SOURCE RADIO MAP was not plotted because they were ALL OUT OF RANGE')
    o0.write('\nNOTE: The SOURCE RADIO MAP was not plotted because they were ALL OUT OF RANGE\n')

print('\n*****************************************')
print('            S  u  m  m  a  r  y')
print('\nWere analyzed %d radio sources' %(nf))
i = nf-rsnp
if (i != 0):
    print(('%d of them were in range: %d with nans + %d without nans' %(i,rsnan,i-rsnan)))
if (rsnid != 0):
    print('The name of %d of them NO were(was) identified' %(rsnid))
if (rsnp != 0):
    print('%d of them were out of range\n' %(rsnp))

o0.write('\n*****************************************\n')
o0.write('            S  u  m  m  a  r  y\n')
o0.write('\nWere analyzed %d radio sources\n' %(nf))
if (i != 0):
    o0.write('%d of them were in range: %d with nans + %d without nans\n' %(nf-rsnp,rsnan,nf-rsnp-rsnan))
if (rsnid != 0):
    o0.write('The name of %d of them NO were(was) identified\n' %(rsnid))
if (rsnp != 0):
    o0.write('%d of them were out of range\n' %(rsnp))
t0 = ctime()
print('\nENDING at %s' %(t0))
o0.write('\nENDING at %s\n' %(t0))
o0.close()
plt.show() 