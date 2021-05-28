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

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

mydata = np.loadtxt('CF_Rx4Ry4Rz4_ellmax4.txt',skiprows=1,unpack=True)
q=mydata[0]
CF_ellmax4=mydata[1]
mydata = np.loadtxt('CF_Rx4Ry4Rz4_ellmax6.txt',skiprows=1,unpack=True)
CF_ellmax6=mydata[1]

xline=[0.0,1000]
yline=[1.0,1.0]
plt.plot(xline,yline,linestyle='--',linewidth=1,color='grey')

plt.plot(q,CF_ellmax4,linestyle='-',linewidth=2,color='b')
#plt.plot(q,CF_ellmax6,linestyle='-',linewidth=2,color='r')



ax.set_xticks(np.arange(0,620,100), minor=False)
ax.set_xticklabels(np.arange(0,620,100), minor=False, family='serif')
ax.set_xticks(np.arange(0,620,50), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
plt.xlim(0.0,230)

ax.set_yticks(np.arange(0,2,0.2), minor=False)
ax.set_yticklabels(np.arange(0,2,0.2), minor=False, family='serif')
ax.set_yticks(np.arange(0,2,0.1), minor=True)
ax.yaxis.set_ticks_position('both')
plt.ylim(0.6,1.2)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$q$ (MeV$/c$)', fontsize=18, weight='normal')
plt.ylabel('$C(q_{\\rm inv})$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('CF_nocoulomb.pdf',format='pdf')
os.system('open -a Preview CF_nocoulomb.pdf')
quit()
