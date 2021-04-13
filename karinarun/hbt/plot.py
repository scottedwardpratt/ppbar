import numpy as np
import matplotlib.pyplot as plt

filename = "CF_1.dat"
mydata = np.genfromtxt(filename, delimiter = " ")[:,:]
x=mydata[:,0]

fig,ax = plt.subplots()

plt.plot(x,mydata[:,1], 'b-', label='weighted')
plt.plot(x,mydata[:,2], 'r-', label='unweighted')

ax.set_title('Weighted CF')
plt.xlabel('q')
plt.ylabel('CF')
#plt.plot(mydata[:,0],mydata[:,2])
#plt.savefig("CF_1.png")
#plt.close()
ax.legend()
plt.show()