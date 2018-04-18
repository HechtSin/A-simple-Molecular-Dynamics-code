import matplotlib.pyplot as plt
import numpy as np

temperature = np.loadtxt('temperature.txt')

plt.plot(temperature)
plt.axhline(y=300,color='r',linestyle='--')
plt.xlim(0,)
plt.ylim(100,500)
plt.xlabel('Timesteps(fs)')
plt.ylabel('Temperature(K)')
#plt.legend(['$\\nu=0.0005$','$\\nu=0.001$','$\\nu=0.002$'])
#plt.title('Leapfrog NVE temperature')
#plt.savefig('leapfrog_nvt_compare.png')
plt.show()
