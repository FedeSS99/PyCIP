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

import numpy as np
from time import ctime,gmtime,localtime,strftime
import os

"""
  The current directory which is where the running program resides,
  is denoted by "." 
 
   At that level is the wipss subdirectory.
 
  And in wipss they are: data, list, spec-data-smooth, spectra-data & 
  spectra-model. Where the data to be analyzed is located and the reports
  and graphs will be stored.
"""


# Conversion factors
deg2rad = np.pi/180.0                  # degrees => radians
rad2deg = 1/deg2rad                    # radians => degrees
# Days of the month
DaysMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
# Spencer solar latitude
SpencerSolLat = lambda g: 0.006918 - 0.399912*np.cos(g)+0.070257*np.sin(g)-0.006758*np.cos(2*g)+0.000907*np.sin(2*g)-0.002697*np.cos(3*g)+0.00148*np.sin(3*g)

class PyCIP:
    def __init__(self, Files:list[str]) -> None:
        # Files to work with
        self.__Files = Files

        # Number of files
        self.__NumFiles = len(Files)

        # Sun coordinates at time of obs
        self.__RaSun = np.zeros(self.__NumFiles)         # solar right ascension (AR) at time of observation
        self.__LaSun = np.zeros(self.__NumFiles)         # solar latitude at the time of observation

        # Source coordinates at time of obs
        self.__RaSouce = np.zeros(self.__NumFiles)       # right ascensions
        self.__DeltaSouce = np.zeros(self.__NumFiles)    # declinations

        # Sun and source coordinates
        self.__Elong = np.zeros(self.__NumFiles)         # elongations
        self.__Arc = np.zeros(self.__NumFiles)           # to plot angles
        self.__HeLat = np.zeros(self.__NumFiles)         # heliocentric latitudes at the time of observation
        self.__DiffRa = np.zeros(self.__NumFiles)        # difference between RA of SOURCE and RA of SUN
        self.__Sides = np.zeros(self.__NumFiles)         # West (W) or East(E) side
        self.__PAngle = np.zeros(self.__NumFiles)        # Position Angle from solar North (same as LASCO)

        self.__Tres = np.zeros((self.__NumFiles, 365))
        self.__TElong = np.zeros((self.__NumFiles, 365))  # elongations along the year
        self.__THeLat = np.zeros((self.__NumFiles, 365))  # heliocentric latitudes along the year

        self.__TSolLats = np.zeros(365)                   # All Solar Latitudes

        # Days per month
        self.__DaysMonth = [31,28,31,30,31,30,31,31,30,31,30,31]

        # Radio source name
        Source = np.empty(10,dtype='U1')
        self.__rsnp = 0                    # Number of Radio Sources NOT Analyzed
        self.__rsnan = 0                   # Number of Radio Sources with NAN
        self.__rsnid = 0                   # Number of Radio Sources NOT Identified

        # Message variable
        self.__message = ""


    def __CreateLogFile(self) -> None:
        self.__path = './wipss/'
        data_dir = self.__path + 'logs/'

        if not os.path.exists(data_dir):
            print(f'Making Logs Dir at {data_dir}...')
            os.makedirs(data_dir)

        LocalTime = strftime("%Y%m%d-%H%M%S",localtime())
        LogName = data_dir + 'log-' + LocalTime + '.txt'
        self.__LogFile = open(LogName, "w")


    def __SolarLatitudesBySpencer(self, index:int, DOY:int):
        g = ((2*np.pi)/365)*(DOY)                 # Gamma diary angle in RADIANS
        lat = SpencerSolLat(g)
        self.__LaSun[index] = lat*rad2deg         # Solar latitude in degrees
        

    
    def StartAnalysis(self) -> None:
        t0 = ctime()

        # Create Log file and print first messages
        self.__CreateLogFile()

        self.__message = '\n     A  U  T  O    -    M  E  X  A  R  T\n'
        print(self.__message)
        self.__LogFile.write(self.__message)
        self.__message = f"Beggining at {t0}"
        print(self.__message)
        self.__LogFile.write(self.__message)

        # Define DATA directory
        self.__datadir = self.__path + "data/" 
        
        self.__message = f"There are {self.__NumFiles} files to process"
        print(self.__message)
        self.__LogFile.write(self.__message)

        # Run over every file
        for nfile in range(self.__NumFiles):
            # Capture current file
            CurrentFile = self.__Files[nfile]

            self.__message = f'\nProcessing file {CurrentFile} - {nfile+1}/{self.__NumFiles}'
            print(self.__message)
            self.__LogFile.write(self.__message)
            self.__message = '(1)  Computing elengation & others'
            print(self.__message)
            self.__LogFile.write(self.__message)

            # Start computation for elongation
            Year, Month, Day = int(CurrentFile[0:4]), int(CurrentFile[4:6]), int(CurrentFile[6:8])
            Date = CurrentFile[0:8]
            RecHour, RecMin, RecSec = int(CurrentFile[9:11]), int(CurrentFile[11:13]), int(CurrentFile[13:15])
            Sign = CurrentFile[15]
            DecDeg, DecMin, DecSec = int(CurrentFile[16:18]), int(CurrentFile[18:20]), int(CurrentFile[20:22])
            if Sign == "-":
                DecDeg *= -1
                DecMin *= -1
                DecSec *= -1

            DOY = 0   # Days accumulated until the begginig of the month. Jan 01 = 0
            for indMonth in range(Month - 1):
                DOY = DaysMonth[indMonth] + DOY

                # Check if leap year
                if (not(Year % 4) and (indMonth == 1) and ((Year%100) or (not(Year%400)))):
                    DOY += 1
            
            DOY += Day

            self.__SolarLatitudesBySpencer(nfile, DOY)
