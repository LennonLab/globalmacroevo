from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import os


# plotting parameters
fs = 10


# modeling parameters
r = 0.0064
Eb = 0.001
K = 10**12



Ss = []
ts = []
S = 1

for t in range(1, 4000):
    ts.append(t)
    S += S*r
    Ss.append(S)


fig = plt.figure()
fig.add_subplot(2, 2, 1)
plt.plot(ts,np.log10(Ss))
plt.text(100, 11, r'$S_{t+1} = S_{t}r$', fontsize=fs-1)
plt.text(100, 8, r'$S_{max}$'+' = '+str('{:.1e}'.format(max(Ss))), fontsize=fs-1)
plt.ylabel(r'$S$', fontsize=fs-1)
plt.xlabel(r'$t$' + ', millions of yrs', fontsize=fs-1)
plt.ylim(0, 13)
plt.xlim(1, 4000)
plt.title('Infinite ' + r'$K$', fontsize=fs)
x = Ss[-1]
print '{:.2e}'.format(x)




Ss = []
ts = []
S = 1

for t in range(1, 4000):
    ts.append(t)
    S += S*r - S*Eb
    Ss.append(S)



fig.add_subplot(2, 2, 2)
plt.plot(ts,np.log10(Ss))
plt.text(100, 11, r'$S_{t+1} = S_{t}r - S_{t}E_{b}$', fontsize=fs-1)
plt.text(100, 8, r'$S_{max}$'+' = '+str('{:.1e}'.format(max(Ss))), fontsize=fs-1)
plt.ylabel(r'$S$', fontsize=fs-1)
plt.xlabel(r'$t$' + ', millions of yrs', fontsize=fs-1)
plt.ylim(0, 13)
plt.xlim(1, 4000)
plt.title('Infinite ' + r'$K$' + ' + BG extinction', fontsize=fs)
x = Ss[-1]
print '{:.2e}'.format(x)



Ss = []
ts = []
S = 1

for t in range(1, 4000):
    ts.append(t)
    S += r*S*((K - S)/K)
    Ss.append(S)


fig.add_subplot(2, 2, 3)
plt.plot(ts,np.log10(Ss))
plt.text(100, 11, r'$S_{t+1} = S_{t} + rS_{t}((K-S_{t})/K)$', fontsize=fs-1)
plt.text(100, 8, r'$S_{max}$'+' = '+str('{:.1e}'.format(max(Ss))), fontsize=fs-1)
plt.ylabel(r'$S$', fontsize=fs-1)
plt.xlabel(r'$t$' + ', millions of yrs', fontsize=fs-1)
plt.ylim(0, 13)
plt.xlim(1, 4000)
plt.title('Finite ' + r'$K$', fontsize=fs)
x = Ss[-1]
print '{:.2e}'.format(x)


Ss = []
ts = []
S = 1

for t in range(1, 4000):
    ts.append(t)
    S += (r-Eb)*S*((K - S)/K)
    Ss.append(S)


fig.add_subplot(2, 2, 4)
plt.plot(ts,np.log10(Ss))
plt.text(100, 11, r'$S_{t+1} = S_{t} + (r-E_{b})S_{t}((K-S_{t})/K)$', fontsize=fs-1)
plt.text(100, 8, r'$S_{max}$'+' = '+str('{:.1e}'.format(max(Ss))), fontsize=fs-1)
plt.ylabel(r'$S$', fontsize=fs-1)
plt.xlabel(r'$t$' + ', millions of yrs', fontsize=fs-1)
plt.ylim(0, 13)
plt.xlim(1, 4000)
plt.title('Finite ' + r'$K$'+ ' + BG extinction', fontsize=fs)
x = Ss[-1]
print '{:.2e}'.format(x)



#plt.show()
plt.suptitle(r'$r$' + ' = ' +str(r) +'  ;  ' +r'$E_{b}$'+' = '+str(Eb)+'  ;  ' +r'$K$'+' = '+str('{:.2e}'.format(K)), fontsize=fs)
plt.subplots_adjust(wspace=0.4, hspace=0.55)


#### Final Format and Save #####################################################
mydir = os.path.expanduser('~/Desktop')
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/ModelFigs.png', dpi=200, bbox_inches = "tight")
plt.close()
