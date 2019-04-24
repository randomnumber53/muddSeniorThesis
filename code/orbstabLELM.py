from itertools import *
from collections import defaultdict

def reduce(state, d):
	return(tuple(map((lambda x: x % d), state)))

def staggerPhase(state, d):
	'''state = (p, c)'''
	(p, c) = state
	return(reduce((p + c, c), d))

def staggerCClass(state, d):
	'''state = (p, c)'''
	(p, c) = state
	return(reduce((p, p + c), d))

def simpleStab(d):
	'''
	Calculates the stabalizer of the element \Psi_{0}^{0} 
	'''
	phase = (1 % d, 0 % d)  # phase             is 1 * p + 0 * c 
	cclass = (0 % d, 1 % d) # correlation class is 0 * p + 1 * c
	initialState = (phase, cclass)
	reachedStates = [initialState]
	currentStates = [initialState]
	while currentStates:
		state = currentStates.pop()
		(currentPhase, currentCClass) = state
		staggeredPhase = (((currentPhase[0] + currentCClass[0]) % d,
						   (currentPhase[1] + currentCClass[1]) % d),
						   currentCClass)
		staggeredCClass = ((currentPhase,
							((currentPhase[0] + currentCClass[0]) % d,
							 (currentPhase[1] + currentCClass[1]) % d)))
		if staggeredPhase not in reachedStates:
			reachedStates.append(staggeredPhase)
			currentStates.append(staggeredPhase)
		if staggeredCClass not in reachedStates:
			reachedStates.append(staggeredCClass)
			currentStates.append(staggeredCClass)
	return(reachedStates)

def multiTest(n):
	for d in range(1, n):
		stab = simpleStab(d)
		stabSize = len(stab)
		print(str(d) + " | " + str(stabSize) + " | " + str(stabSize * d**2))

def groupSize(d):
	return(len(simpleStab(d)*d*d))

def stabalizer(subbasis, d):
	'''
	subbasis is a set of tuples, which each tuple (p, c) represents the
	state \Psi_{c}^{p}
	'''
	subbasis = set(subbasis)
	group = []
	out = []
	# First, we construct the group (this need only be done once per d)
	for ss in simpleStab(d):
		for i in range(d):
			for j in range(d):
				group.append((ss[0] + (i,), ss[1] + (j,)))
	for g in group:
		transform = (lambda q: ((g[0][0] * q[0] + g[0][1] * q[1] + g[0][2]) % d,
								(g[1][0] * q[0] + g[1][1] * q[1] + g[1][2]) % d))
		newSubbasis = set(map(transform, subbasis))
		if newSubbasis == subbasis:
			out.append(g)

	return(sorted(out))

def bellStates(d):
	return(list(product(range(d), repeat=2)))

def subbases(d, k):
	return(list(combinations(bellStates(d), k)))

## This creates a histogram of the sizes of orbits, which lets us calculate the
## number of orbits of each size.

# a = [len(stabalizer(b, 4)) for b in subbases(4, 7)]
# d = defaultdict(int)
#for n in a: d[n] += 1 

