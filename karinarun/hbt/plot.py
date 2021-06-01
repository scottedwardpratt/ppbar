import numpy as np
import matplotlib.pyplot as plt

filename = "exercise.dat"
mydata = np.genfromtxt(filename, delimiter = " ")[:,:]
ix=mydata[:,0]

#fig,ax = plt.subplots()

plt.plot(ix,mydata[:,1], 'o')
#plt.plot(x,mydata[:,2], 'r-', label='unweighted')

#ax.set_title('Boltzmann')
plt.xlabel('ix')
plt.ylabel('dn/dx')
#plt.plot(mydata[:,0],mydata[:,2])
#plt.savefig("CF_1.png")
#plt.close()
#ax.legend()
plt.show()