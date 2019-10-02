from __future__ import division
import  matplotlib.pyplot as plt
from numpy.random import uniform
import numpy as np
import os
import sys

sz = 1
def assigncolor(xs):
    clrs = []

    for x in xs:
        c = 0
        if x <= 6: c = 'crimson'
        elif x <= 7: c = 'r'
        elif x <= 8: c = 'orangered'

        elif x <= 9: c = 'darkorange'
        elif x <= 10: c = 'orange'
        elif x <= 11: c = 'gold'

        elif x <= 12: c = 'yellow'
        elif x <= 13: c = 'greenyellow'
        elif x <= 14: c = 'springgreen'

        elif x <= 15: c = 'Green'
        elif x <= 16: c = 'skyblue'
        elif x <= 17: c = 'DodgerBlue'

        elif x <= 18: c = 'b'
        elif x <= 19: c = 'Plum'
        elif x <= 20: c = 'violet'
        elif x <= 21: c = 'purple'
        else: c = 'k'

        clrs.append(c)
    return clrs



# plotting parameters
fs = 10

rs = []
Ebs = []
Ks = []
Ss = []

for i in range(40000):
    # modeling parameters
    r = uniform(0.002, 0.03)
    rs.append(r)

    Eb = 10**uniform(-3, 0) # epsilon
    Ebs.append(np.log10(Eb))

    Eb = Eb*r # convert epsilon to mu
    K = uniform(0, 1)
    Ks.append(np.log10(K))

    S = 1
    for t in range(1, 4000):
        S += (r-Eb)*S*((K - S)/K)
    Ss.append(np.log10(S))


clrs = assigncolor(Ss)

fig = plt.figure()
fig.add_subplot(3, 3, 1)
plt.scatter(rs, Ks, color=clrs, s=sz, linewidths=0.0, edgecolor=None, alpha=0.5)
plt.xlabel(r'$r$', fontsize=fs)
plt.ylabel(r'$K$', fontsize=fs)
plt.xlim(min(rs), max(rs))
plt.ylim(min(Ks), max(Ks))
plt.tick_params(axis='both', labelsize=fs-4)



fig.add_subplot(3, 3, 2)
sz = 16

plt.scatter([1], [1], color='crimson', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{6}$')
plt.scatter([1], [1], color='r', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{7}$')
plt.scatter([1], [1], color='orangered', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{8}$')
plt.scatter([1], [1], color='darkorange', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{9}$')
plt.scatter([1], [1], color='orange', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{10}$')
plt.scatter([1], [1], color='gold', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{11}$')
plt.scatter([1], [1], color='yellow', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{12}$')
plt.scatter([1], [1], color='greenyellow', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{13}$')
plt.scatter([1], [1], color='springgreen', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{14}$')
plt.scatter([1], [1], color='Green', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{15}$')
plt.scatter([1], [1], color='skyblue', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{16}$')
plt.scatter([1], [1], color='b', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{17}$')
plt.scatter([1], [1], color='DodgerBlue', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{18}$')
plt.scatter([1], [1], color='Plum', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{19}$')
plt.scatter([1], [1], color='violet', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{20}$')
plt.scatter([1], [1], color='purple', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{21}$')
plt.scatter([1], [1], color='k', s=sz, linewidths=0.0, edgecolor=None, label=r'$S < 10^{21}$')

plt.ylim(2,3)
plt.xlim(2,3)
plt.xticks([])
plt.yticks([])
plt.box(on=None)
plt.legend(loc=2, fontsize=fs, frameon=False)
sz = 1

fig.add_subplot(3, 3, 4)
plt.scatter(rs, Ebs, color=clrs, s=sz, linewidths=0.0, edgecolor=None, alpha=0.5)
plt.xlabel(r'$r$', fontsize=fs)
plt.ylabel(r'$E_{b}$', fontsize=fs)
plt.xlim(min(rs), max(rs))
plt.ylim(min(Ebs), max(Ebs))
plt.tick_params(axis='both', labelsize=fs-4)



fig.add_subplot(3, 3, 7)
plt.scatter(Ebs, Ks, color=clrs, s=sz, linewidths=0.0, edgecolor=None, alpha=0.5)
plt.xlabel(r'$E_{b}$', fontsize=fs)
plt.ylabel(r'$K$', fontsize=fs)
plt.xlim(min(Ebs), max(Ebs))
plt.ylim(min(Ks), max(Ks))
plt.tick_params(axis='both', labelsize=fs-4)



#plt.show()
plt.subplots_adjust(wspace=0.4, hspace=0.4)


#### Final Format and Save #####################################################
mydir = os.path.expanduser('~/Desktop')
plt.subplots_adjust(wspace=0.4, hspace=0.5)
plt.savefig(mydir + '/ModelFigs.png', dpi=200, bbox_inches = "tight")
plt.close()
