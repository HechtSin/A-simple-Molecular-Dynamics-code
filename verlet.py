import numpy as np
from LJ import potential,force
from initialize import total_num, position, velocity,lc,mass,Temp,box
import matplotlib.pyplot as plt
from scipy.constants import k
import random

###     Main FILE      ###
### Verlet Algorithm ###

### INPUTS by users ##########################################################

## step ##
h = 1e-15 ## 1 femtosecond
## how many steps we wanna run
total_step = 500
## total time
tmax = h*total_step

## A variable to describe coupling between system and heat bath (only for NVT)
nu = 0.001

##################################################################
## Initial temperature
temperature = 0
## Initial kinetic energy
KE0 = 0
## Initial potential energy
PE0 = 0

## Get force at t = 0 
Force = np.zeros_like(position)
for i in range(total_num):
    for j in range(total_num):
        if j != i:
            Force[i] += force(position[i],position[j])

## Calculate initial energies
for i in range(total_num):
    temperature += 1/3*mass*np.linalg.norm(velocity[i])**2/k/total_num
    KE0 += 0.5*mass*np.linalg.norm(velocity[i])**2
    for j in range(total_num):
        if j!=i:
            PE0 += potential(position[i],position[j])/2
TE0 = KE0 + PE0

print ('Verlet algorithm-NVE system')

print ('KE0={},PE0={},TE0={}'.format(KE0,PE0,TE0))

print ('Initial temperature = {}K'.format(temperature))


####################################################################


## Define variables to store information at half-integer steps
position_half = np.zeros_like(position)
velocity_half = np.zeros_like(position)
Force = np.zeros_like(position)
Force_half = np.zeros_like(position)

## Get velocity at 0.5h step
velocity_half = velocity + 0.5*h*Force/mass

## Run MD verlet

## lists to append energy and temperature
energy = []
all_temp = []


### main  loop
t = h
while t<tmax:

    ## Print to let user know where we are
    if round(t/h)%100 ==0:
        print ('Step {}'.format(round(t/h)))

    ## position at t+h
    position = position + h*velocity_half

    ## Force at t+h
    Force = np.zeros_like(position)
    for i in range(total_num):
        for j in range(total_num):
            if j != i:
                Force[i] += force(position[i],position[j])

    ## a temporary variable
    temporary = h*Force/mass
    
    ## velocity at t+h
    velocity = velocity_half + 0.5*temporary

    ## velocity at t+1.5h
    velocity_half = velocity_half + temporary

    ## calc energy and temperature
    KE = 0
    PE = 0
    TE = 0
    temperature = 0
    for i in range(total_num):
        for j in range(total_num):
            if j != i:
                PE += potential(position[i],position[j])

        KE += 0.5*mass*np.linalg.norm(velocity[i])**2
        temperature += 1/3*mass*np.linalg.norm(velocity[i])**2/k/total_num
  
    TE = KE +PE 

    ## Andersen thermostat
    ## For NVE, commented out this part 
    #for i in range(total_num):
    #    a = random.uniform(0,1)
    #    if a < nu:
    #        velocity[i][0] = np.asscalar(np.random.normal(loc=0.0,scale=1,size=1)*np.sqrt(k*Temp/mass))
    #        velocity[i][1] = np.asscalar(np.random.normal(loc=0.0,scale=1,size=1)*np.sqrt(k*Temp/mass))
    #        velocity[i][2] = np.asscalar(np.random.normal(loc=0.0,scale=1,size=1)*np.sqrt(k*Temp/mass))

    energy.append(TE)
    all_temp.append(temperature)
   
    ## If someting is wrong, print out 
    if TE >0: ## energy is larger than 0, somethign wrong
       print (position)
       print (velocity)
       break

    t += h

### save files for further plotting
energy = np.array(energy)
all_temp = np.array(all_temp)

np.savetxt('energy.txt',energy)
np.savetxt('temperature.txt',all_temp)
