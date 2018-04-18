import matplotlib.pyplot as plt
import numpy as np

energy = np.loadtxt('energy.txt')

plt.plot(energy)
plt.xlim(0,)
plt.xlabel('Timesteps(fs)')
plt.ylabel('Total energy (J)')
#plt.title('Verlet NVE energy')
#plt.title('Leapfrog NVT energy')
plt.gcf().subplots_adjust(left=0.15)
#plt.savefig('leapfrog_nve_energy.png')
#plt.savefig('leapfrog_nvt_energy.png')
plt.show()
