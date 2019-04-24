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

def stab(d):
	phase = (1 % d, 0 % d)
	cclass = (0 % d, 1 % d)
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
		stab = test(d)
		stabSize = len(stab)
		print(str(d) + " | " + str(stabSize) + " (" + str(d**4) + ", " + str(d**3-d) + ")" + " | " + str(stabSize * d**2))
