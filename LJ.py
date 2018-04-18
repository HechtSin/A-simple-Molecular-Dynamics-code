import numpy as np
from scipy.constants import k
from initialize import lc,mass,box
import matplotlib.pyplot as plt


###         Potential function FILE            ## 
###        Lennard-Johnes potential            ##
### Describe interaction between any two atoms ##


## Potential function
## Return potential energy Vij
def potential(ri,rj):
  
    ## Get atomic position 
    xi = ri[0]
    yi = ri[1]
    zi = ri[2]

    xj = rj[0]
    yj = rj[1]
    zj = rj[2]
  
    xij = ri[0]-rj[0]
    yij = ri[1]-rj[1]
    zij = ri[2]-rj[2]

    ## Lennard-Jones Potential equation
    ## V_ij(r) = 4*epsilon*((sigma/r)^12-(sigma/r)^6)
    ## F_ij(r) = -12/|r^14|*r_vec+6/|r^8|*r_vec

    ## A parameter that is close to equilibrium position
    sigma = lc/np.sqrt(2)/1.122
    ## Potential depth
    epsilon = 6767.4*k

    ## Periodic boundary conditions
    if xij > box/2:
        xij -= box
    elif xij < -box:
        xij += box

    if yij > box/2:
        yij -= box
    elif yij < -box/2:
        yij += box
    
    if zij > box/2:
        zij -= box
    elif zij < -box/2:
        zij += box

    ## Distance between i,j
    rij = np.sqrt(xij**2+yij**2+zij**2)

    ## Calculate potential
    Vij = 4*epsilon*((sigma/rij)**12-(sigma/rij)**6)
    
    return Vij

## Test
## To plot V vs r
## Should see a potential function similar in slides
#sigma = lc/np.sqrt(2)/1.122
#epsilon = 6767.4*k
#r = np.linspace(0.5*lc,2*lc,100)
#pot = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
#plt.plot(r,pot)
#plt.axhline(y=0.0)
#plt.axvline(x=lc/np.sqrt(2))
#plt.xlim(0.5*lc,1.5*lc)
#plt.ylim(-1e-19,1e-19)
#plt.show()




## Calculate Forces
## Return force between i,j

def force(ri,rj):
   
    ## Get atomic position 
    xi = ri[0]
    yi = ri[1]
    zi = ri[2]

    xj = rj[0]
    yj = rj[1]
    zj = rj[2]
  
    xij = ri[0]-rj[0]
    yij = ri[1]-rj[1]
    zij = ri[2]-rj[2]

    ## Lennard-Jones Potential equation
    ## V_ij(r) = 4*epsilon*((sigma/r)^12-(sigma/r)^6)
    ## F_ij(r) = -12/|r^14|*r_vec+6/|r^8|*r_vec

    ## A paramete close to equilibrium position
    sigma = lc/np.sqrt(2)/1.122
    ## Potential depth
    epsilon = 6767.4*k
    
    ## Periodic boundary conditions 
    if xij > box/2:
        xij -= box
    elif xij < -box/2:
        xij += box

    if yij > box/2:
        yij -= box
    elif yij < -box/2:
        yij += box
    
    if zij > box/2:
        zij -= box
    elif zij < -box/2:
        zij += box
   
    ## Distance between i,j 
    rij = np.sqrt(xij**2+yij**2+zij**2)

    ## Negative derivative of potential
    factor = 4*epsilon*(sigma**12*(12)/(rij**13)+sigma**6*(-6)/(rij**7))

    ## Force component for x,y,z direction
    Fijx = factor*xij/rij
    Fijy = factor*yij/rij
    Fijz = factor*zij/rij

    return np.array([Fijx,Fijy,Fijz])

## Test ##
## To plot force vs r
## Should find force=0,when r is at lowest point of potential

#sigma = lc/np.sqrt(2)/1.122
#epsilon = 6767.4*k
#r = np.linspace(0.5*lc,2*lc,100)
#f = 4*epsilon*(sigma**12*(12)/(r**13)+sigma**6*(-6)/(r**7))
#plt.plot(r,f,'r')
#plt.xlim(0.5*lc,1.5*lc)
#plt.ylim(-1e-8,1e-8)
#plt.show()
