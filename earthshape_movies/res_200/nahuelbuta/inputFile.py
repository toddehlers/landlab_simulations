"""
This file contains the input-parameters for the landlab driver. 

This does NOT use the landlab inherent load_params method but declares
variables directly in python and is loaded in the driver file with

	from inputFile.py import *

This was done because we do not use the standard-inputs for the
fluvial and hillslope routines but do the processing for
vegetation influence directly in the driver.

Usage:

	-Parameternames must be equal to the ones declared in the 
	driver file (better not change at all.)
	-Comments can be written in a standart python way
	-Any declaration of the same variable in the driver-file
	will overwrite the imported values so code-with-caution.


Created by: Manuel Schmid, May 28th, 2018
"""

# Nahuelbuta

#Model Grid Parameters
ncols = 200	#number of columns
nrows = 200  #number of rows
dx    = 200 #spacing between nodes

#Model Runtime Parameters
totalT = 10000000.0  #total model runtime
ssT    = 10000000.0  #spin-up time before sin-modulation, set to same value as totalT for steady-state-simulations
sfT    = 10000000.0  #spin-up time before step-change-modulation, set to same value as totalT for steady-state-simulations
dt     =    10000.0
#Uplift
upliftRate = 0.00001 #m/yr, Topographic uplift rate

#Surface Processes
#Linear Diffusion:
linDiffBase = 3e-1 #m2/yr, base linear diffusivity for bare-bedrock
alphaDiff   = 0.3  #Scaling factor for vegetation-influence (see Instabulluoglu and Bras 2005)

#Fluvial Erosion:
ksp = 2e-7 #base fluvial erodibility for bare-bedrock
msp = 0.5  #m factor from SPL
nsp = 1    #n factor from SPL
thresholdSP = 0.0 #threshold erosion-factor from SPL
critArea    = 1e6 #L^2, Minimum Area which the steepness-calculator assumes for channel formation.
aqDens      = 1000 #Kg/m^3, density of water
grav        = 9.81 #m/s^2, acceleration of gravity
nSoil       = 0.01 #Mannings number for bare soil
nVRef       = 0.6  #Mannings number for reference vegetation
vRef        = 1    #1 = 100%, reference vegetation-cover for fully vegetated conditions
w           = 1    #Scaling factor for vegetation-influence (see Istanbulluoglu and Bras 2005)

#Climate Parameters
baseRainfall = 140.0 #m/dt, base steady-state rainfall-mean over the dt-timespan
rfA          = 0  #m, rainfall-step-change if used
maxRain      = 60 #m/dt, IF you use the rainfall oscillation, -> upper amplitude
lowRain      = 15 #s.o -> lower amplitude

#Vegetation Cover
vp = 0.9 #initial vegetation cover, 1 = 100%
sinAmp = 0.3 #vegetation cover amplitude for oscillation
sinPeriod = 1e6 #yrs, period of sin-modification

#output
outInt = 100000.0 #yrs, model-time-interval in which output is created

