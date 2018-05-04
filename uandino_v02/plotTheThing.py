#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def normalize(data):
	tot = data[:,1]+data[:,2]+data[:,3]
	data[:,1]/=tot
	data[:,2]/=tot
	data[:,3]/=tot


probData=np.loadtxt('probsTest.csv', delimiter=',', dtype=float)
probabilities, ax = plt.subplots(2, 2, figsize=(10, 5))
#normalize(probData)
ax[0, 0].plot(probData[:,0]*1e-6, probData[:,1])
ax[0, 1].plot(probData[:,0]*1e-6, probData[:,2])
ax[1, 0].plot(probData[:,0]*1e-6, probData[:,3])
ax[1, 1].plot(probData[:,0]*1e-6, probData[:,1]+ probData[:,2]+ probData[:,3] )
ax[0, 0].set_ylabel('$P_{e e}$', fontsize=15)
ax[0, 1].set_ylabel('$P_{\mu e}$', fontsize=15)
ax[1, 0].set_ylabel('$P_{\\tau e}$', fontsize=15)
ax[1, 1].set_ylabel('$P_{\\tau e}+P_{\mu e}+ P_{e e}$', fontsize=15)
ax[1, 1].set_ylim(1-1, 1+1)
#ax[0,0].errorbar(np.array([0.3e6, 1e6, 1.1e7])*1e-6,[0.82, 0.26, 0.27], yerr=[0.14, 0.16, 0.06], marker='o')
for i in range(2):
    for k in range(2):
        ax[i,k].set_xscale('log')
        ax[i,k].set_xlabel('$E_{\\nu}$(MeV)', fontsize=15)
        ax[i,k].set_xlim(1e3*1e-6, 1e13*1e-6)
#ax[1,1].set_ylim(1-0.00001, 1+0.00001)
ax[0,1].set_ylim(0,1)
ax[1,0].set_ylim(0,1)
ax[0,0].set_ylim(0,1)

plt.tight_layout()
plt.gcf()
plt.savefig('probPlot.png', dpi=300)

newfig = plt.figure()
pathData = np.loadtxt('potentialTest.csv', delimiter=',', dtype=float)
plt.plot(pathData[:,0], pathData[:,1])
plt.gcf()
#plt.gca().set_yscale('log')
plt.savefig('potentialPlot.png', dpi=300)
