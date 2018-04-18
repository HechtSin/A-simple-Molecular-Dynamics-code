import numpy as np
from scipy.constants import pi,k
import matplotlib.pyplot as plt
from scipy.stats import maxwell


####              Initializing FILE             ##########
###          Read the initial position file     ##########
### Generate initial velocity based on random distribution

## lists to store position variables 
positionx = []
positiony = []
positionz = []

## lattice constants
lc = 3.6150500774*1e-10 # Angstrom

## box size
## 2*2*2 supercell
box = 2*lc

## read initial position file
with open('ssposcar.txt','r') as supercell:
     lines = supercell.readlines()[8:]
     total_num = int(len(lines)) # Total atom number
     for i in lines:
         positionx.append(float(i.split()[0])*box)
         positiony.append(float(i.split()[1])*box)
         positionz.append(float(i.split()[2])*box)
supercell.close

position = np.stack((positionx,positiony,positionz)).T

### Introduce initial velocities
### Use numpy to generate random numbers based on Maxwell-Boltzmann distribution

mass = 1.0552061e-25 ## copper atom
Temp = 300 ## Kelvin

velocityx = np.random.normal(loc=0.0,scale=1,size=total_num)*np.sqrt(k*Temp/mass)
velocityy = np.random.normal(loc=0.0,scale=1,size=total_num)*np.sqrt(k*Temp/mass)
velocityz = np.random.normal(loc=0.0,scale=1,size=total_num)*np.sqrt(k*Temp/mass)
velocity = np.stack((velocityx,velocityy,velocityz)).T

### Calculate net momentum
p = np.zeros(3)
for i in range(total_num):
    p += mass*velocity[i]

### Subtract net momentum
### Whole system should have no net momentum
for i in range(total_num):
    velocity[i] -= p/(total_num*mass)

### Check if net momtentum is now zero
#p = np.zeros(3)
#for i in range(total_num):
#    p += mass*velocity[i]
#print (p)
