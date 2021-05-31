import numpy as np
import matplotlib.pyplot as plt

filename = "test.dat"
mydata = np.genfromtxt(filename, delimiter = " ")[:,:]
ix=mydata[:,0]

fig,ax = plt.subplots()

plt.plot(ix,mydata[:,1])
#plt.plot(x,mydata[:,2], 'r-', label='unweighted')

ax.set_title('Boltzmann')
plt.xlabel('ix')
plt.ylabel('dx/dx')
#plt.plot(mydata[:,0],mydata[:,2])
#plt.savefig("CF_1.png")
#plt.close()
ax.legend()
plt.show()