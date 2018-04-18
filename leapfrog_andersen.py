import numpy as np
from LJ import potential,force
from initialize import total_num, position, velocity,lc,mass,Temp
import matplotlib.pyplot as plt
from scipy.constants import k
import random

###     Main FILE      ###
### Leapfrog Algorithm ###


### INPUTS by users ##########################################################

## step ##
h = 1e-15 ## 1 femtosecond
## how many steps we wanna run
total_step = 500
## total time
tmax = h*total_step

## A variable to describe coupling between system and heat bath (only for NVT)
nu = 0.001



##############################################################################
## Define derivative funciton
## Input position and velocity
## Output velocity and accelerate

def f(r,t):
   ## Position and velocity at time t
   position = r[0]
   velocity = r[1] 
   
   fposition = velocity
   
   fvelocity = np.zeros_like(fposition)

   for i in range(total_num):
       for j in range(total_num):
           if j != i:
               fvelocity[i] += force(position[i],position[j])/mass

   ## Return velocity and acceleration
   return np.stack((fposition,fvelocity))

## initial position is at t=0
## initial velocity is at t=0

## Initial temperature
temperature = 0
## Initial kinetic energy
KE0 = 0
## Initial potential energy
PE0 = 0

## Calcualte initial energies
for i in range(total_num):
    temperature += 1/3*mass*np.linalg.norm(velocity[i])**2/k/total_num
    KE0 += 0.5*mass*np.linalg.norm(velocity[i])**2
    for j in range(total_num):
        if j!=i:
            PE0 += potential(position[i],position[j])/2
TE0 = KE0 + PE0

print ('Leapfrog algorithm')

print ('KE0={},PE0={},TE0={}'.format(KE0,PE0,TE0))

print ('Initial temperature = {}'.format(temperature))


###########################################################################

## Run MD leapfrog ##

## starting time
t = 0

## Used to append total energy and temperature
energy = []
total_temp = []

### Define variables to store information at half-integer steps
position_half = np.zeros_like(position)
velocity_half = np.zeros_like(position)
Force = np.zeros_like(position)
Force_half = np.zeros_like(position)

## initial force
Force = np.zeros_like(position)
for i in range(total_num):
    for j in range(total_num):
        if j != i:
            Force[i] += force(position[i],position[j])

r = np.stack((position,velocity))
### Now, all the variables at t=0 are determined


## first step
## t = 0.5h
for i in range(total_num):
    position_half[i] = position[i] + 0.5*h*velocity[i]
    velocity_half[i] = velocity[i] + 0.5*h*Force[i]/mass

Force_half = np.zeros_like(position)
for i in range(total_num):
    for j in range(total_num):
        if j!= i:
            Force_half[i] += force(position_half[i],position_half[j])

r_half = np.stack((position_half,velocity_half))

### Now, all the variables at t=1/2h are dtermined


#### main loop ####

t = h
while t<tmax:
    ## Print to make user know where we are
    if round(t/h)%100 ==0:
        print ('Step {}'.format(round(t/h)))

    ## variables at t+h
    r += h*f(r_half,t)
    ## variables at t+3/2h
    r_half += h*f(r,t)

    ## calc energy and temperature at t+h
    KE = 0
    PE = 0
    TE = 0
    temperature = 0
    for i in range(total_num):
        for j in range(total_num):
            if j != i:
                PE += potential(r[0][i],r[0][j])

        KE += 0.5*mass*np.linalg.norm(r[1][i])**2
        temperature += 1/3*mass*np.linalg.norm(r[1][i])**2/k/total_num
  
    TE = KE +PE 

    energy.append(TE)
    total_temp.append(temperature) 
   
    ## Andersen thermostat
    ## If use NVE, comment out this part
    ## After each step, there is probability to collide
    ## If collides, the velocity is re-distributed by Maxwell-Boltzmann
    for i in range(total_num):
        a = random.uniform(0,1)
        if a < nu:
            r[1][i][0] = np.asscalar(np.random.normal(loc=0.0,scale=1,size=1)*np.sqrt(k*Temp/mass))
            r[1][i][1] = np.asscalar(np.random.normal(loc=0.0,scale=1,size=1)*np.sqrt(k*Temp/mass))
            r[1][i][2] = np.asscalar(np.random.normal(loc=0.0,scale=1,size=1)*np.sqrt(k*Temp/mass))


    ## If something is wrong
    ## break to check
    if TE >0: ## somthing is wrong, if energy>0
       print (position)
       print (velocity)
       break

    t += h

## Save files for further plotting
energy = np.array(energy)
total_temp = np.array(total_temp)
np.savetxt('temperature.txt',total_temp)
np.savetxt('energy.txt',energy)
