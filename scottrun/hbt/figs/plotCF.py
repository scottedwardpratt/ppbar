import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))
pi=np.pi
import sys

#print(sys.argv[1])

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(8,6))
fig = plt.figure(1)
ax = fig.add_axes([0.125,0.125,0.84,0.84])

#filename='../results_'+sys.argv[1]+'/CFs/'+sys.argv[2]+'.txt'
filename= "../results_qcut/CFs/average.txt"
print('filename=',filename)

mydata = np.loadtxt(filename,skiprows=2,unpack=True)
plt.plot(mydata[0],mydata[1],linestyle='-',linewidth=2,color='k')
plt.plot(mydata[0],mydata[2],linestyle='-',linewidth=2,color='r')
plt.plot(mydata[0],mydata[3],linestyle='-',linewidth=2,color='g')
plt.plot(mydata[0],mydata[4],linestyle='-',linewidth=2,color='b')

ax.tick_params(axis='both', which='major', labelsize=18)

ax.set_xticks(np.arange(0,120,20), minor=False)
ax.set_xticklabels(np.arange(0,120,20), minor=False, family='serif')
ax.set_xticks(np.arange(0,120,5), minor=True)
plt.xlim(0,100)

ax.set_yticks(np.arange(0,2.0,0.5), minor=False)
ax.set_yticklabels(np.arange(0,2.0,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(0,2.0,0.1), minor=True)
plt.ylim(0,1.0)

plt.xlabel('$q$ (MeV/c)', labelpad=0, fontsize=24, weight='normal')
plt.ylabel('$C(q)$',fontsize=24)
plt.savefig('CF.pdf',format='pdf')
os.system('open -a Preview CF.pdf')

#plt.show()
quit()
