# Simulate variable rates among clades
# Ford Fishman

from math import *


S_total = [1.0] # initial total richness
lam = 0.015 # speciation rate
epsilon = [float(x)/20 for x in range(2, 19)] # extinction/speciation ranging from 0.1 - 0.9
species = {ep:[0.0] for ep in epsilon}
species[0.5] = [1.0]
mu = [ep * lam for ep in epsilon]
r = [lam - m for m in mu]


def timestep(ep, clades, i):

	if i == 0:
		clades[ep].append(0.97 * clades[ep][-1] * (r[i]) + 
			0.02 * clades[epsilon[i+1]][-1] * (r[i+1]) + 
			0.01 * clades[epsilon[i+2]][-1] * (r[i+2]) +
			clades[ep][-1])
	elif i == 1:
		clades[ep].append(0.95 * clades[ep][-1] * (r[i]) + 
			0.02 * clades[epsilon[i-1]][-1] * (r[i-1]) + 
			0.02 * clades[epsilon[i+1]][-1] * (r[i+1]) + 
			0.01 * clades[epsilon[i+2]][-1] * (r[i+2]) +
			clades[ep][-1])
	elif epsilon[i] == epsilon[-2]:
		clades[ep].append(0.95 * clades[ep][-1] * (r[i]) + 
			0.02 * clades[epsilon[i-1]][-1] * (r[i-1]) + 
			0.02 * clades[epsilon[i+1]][-1] * (r[i+1]) + 
			0.01 * clades[epsilon[i-2]][-1] * (r[i-2]) +
			clades[ep][-1]) 
	elif epsilon[i] == epsilon[-1]:	
		clades[ep].append(0.97 * clades[ep][-1] * (r[i]) + 
			0.02 * clades[epsilon[i-1]][-1] * (r[i-1]) + 
			0.01 * clades[epsilon[i-2]][-1] * (r[i-2]) +
			clades[ep][-1]) 
	else:
		clades[ep].append(0.94 * clades[ep][-1] * (r[i]) + 
			0.02 * clades[epsilon[i-1]][-1] * (r[i-1]) + 
			0.02 * clades[epsilon[i+1]][-1] * (r[i+1]) + 
			0.01 * clades[epsilon[i-2]][-1] * (r[i-2]) + 
			0.01 * clades[epsilon[i+2]][-1] * (r[i+2]) +
			clades[ep][-1])
	return(clades)

for i in range(1, 4000):
	S = []
	for ep in species:
		species = timestep(ep = ep, clades = species, i = epsilon.index(ep))
		S.append(species[ep][-1])
	S_total.append(sum(S))

	# species[0.5].append(0.94 * species[0.5][-1] * (r[1]) + 
	# 	0.03 * species[0.1][-1] * (r[0]) + 
	# 	0.03 * species[0.9][-1] * (r[2])) 
	# species[0.1].append(0.03 * species[0.5][-1] * (r[1]) + 
	# 	0.97 * species[0.1][-1] * (r[0]) + 
	# 	0 * species[0.9][-1] * (r[2])) 
	# species[0.9].append(0.03 * species[0.5][-1] * (r[1]) + 
	# 	0 * species[0.1][-1] * (r[0]) + 
	# 	0.97 * species[0.9][-1] * (r[2]))
	# S_total.append(species[0.1][-1] + species[0.5][-1] + species[0.9][-1])

for ep in species:
	print(str(ep) + ": "+ str(species[ep][-1]))
print("Total: " + str(S_total[-1]))