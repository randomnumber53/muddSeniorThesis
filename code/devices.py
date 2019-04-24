from itertools import *
from math import *
from cmath import *

cdev = [(0, 0, 0, 0, 1, 0), (0, 0, 1, 0, 0, 0), (0, 0, 1, 0, 1, 0),
        (0, 0, 1, 1, 0, 0), (0, 0, 1, 1, 1, 0), (0, 1, 0, 1, 0, 0)]

d = 4
omega = rect(1, 2 * pi / d)

def innerProduct(left, right):
    return sum(p.conjugate()*q for p, q in zip(left, right))

def expandBellState(compState):
    c, p = compState
    state = [0]*(d**2)
    for i in range(d):
        state[(i*d + c) % d**2] = omega**(p * i)
    return tuple(state)

def buildDevices(k):
    return list(combinations(signatures, k))

def getSignature(mode1, mode2):
    """
    (|00>,  |01>,  |02>,
     |10>,  |11>,  |12>,
     |20>,  |21>,  |22>)
    """
    dim = len(mode1)//2
    result = [0]*(dim**2)
    for i in range(dim):
        for j in range(dim):
            result[i + dim*j] = mode1[i]*mode2[j+dim] + mode2[i]*mode1[j+dim]
    return(tuple(result))

def testDevice(dev, states):
    """
    Returns True if dev can distinguish the states and false otherwise.
    """
    spansStates = [False]*len(states) 
    noReps = True

    for i in range(1, len(dev)):
        currentMode = modes[(dev[0], dev[i])]
        modeUsed = False
        rep = False
        for l in range(len(states)):
            if polar(innerProduct(currentMode, states[l]))[0] > 0.0001:
                modeUsed = True
                spansStates[l] = True
                if rep == True:
                    return False
                rep = True
        if not modeUsed:
            return False

    return sum(spansStates) == len(states)
        
def buildModes():
    modes = {}
    divisor = len(signatures) // 1000
    for i in range(len(signatures)):
        for k in range(i+1, len(signatures)):
            modes[(signatures[i], signatures[k])] = getSignature(signatures[i],
                                                            signatures[k])
        if i % divisor == 0:
            print(str(100 * round(i / len(signatures), 4)) + "% of modes built")
    return(modes)

def testAllDevices(k, states):
    devs = buildDevices(k)
    divisor = len(devs) // 1000
    for i in range(len(devs)):
        dev = devs[i]
        if testDevice(dev, states):
            print(states)
            print(dev)
            return(True)
        if i % divisor == 0:
            print(str(100 * round(i / len(devs), 4)) + "% of devices checked")
    return(False)

def analyzeDevice(dev, states):
    statesToSigs = {state: [] for state in states}
    print("The device includes the modes:")
    for i in range(len(dev)):
        print(dev[i])
    print("They omit the detection signatures:")
    for i in range(len(dev)):
        for k in range(i+1, len(dev)):
            currentMode = modes[(dev[i], dev[k])]
            print(currentMode)
            for state in states:
                if polar(innerProduct(currentMode, state))[0] > 0.0001:
                    statesToSigs[state].append((i, k))
    print("Which clicks correspond to which states?")
    for i in range(len(states)):
        print("state " + str(i) + ": " + str(statesToSigs[states[i]]))



# |0,L> + |1,L> + |2,L> + |0,R> + |1,R> + |2,R>
modes1 = [
(1, 0, 1, 0, 1, 0, 1, 0),
(1, 0, 1j, 0, -1, 0, -1j, 0),
(1, 0, -1, 0, 1, 0, -1, 0),
(1, 0, -1j, 0, -1, 0, -1j, 0),
(0, 1, 0, 1, 0, 1, 0, 1),
(0, 1, 0, 1j, 0, -1, 0, -1j),
(0, 1, 0, -1, 0, 1, 0, -1),
(0, 1, 0, -1j, 0, -1, 0, -1j)
]

modes2 = [
(1, 0, 0, 0, 1, 0, 0, 0),
(0, 1, 0, 0, 0, 1j, 0, 0),
(0, 0, 1, 0, 0, 0, -1, 0),
(0, 0, 0, 1, 0, 0, 0, -1j),
(1, 0, 0, 0, 1j, 0, 0, 0),
(0, 1, 0, 0, 0, -1, 0, 0),
(0, 0, 1, 0, 0, 0, -1j, 0),
(0, 0, 0, 1, 0, 0, 0, 1)
]

modes = modes2
# (c, p)
compBellBasis = list(product(range(d), repeat = 2))

bellBasis = list(map(expandBellState, compBellBasis))

testBasis1 = [bellBasis[0], bellBasis[4], bellBasis[5], bellBasis[6], bellBasis[14]]

allBases = list(combinations(bellBasis, 5))

def test(testBasis = testBasis1, p=False):
    largest = 0
    for i in range(len(modes)):
        for j in range(i, len(modes)):
            statesPresent = []
            sig = getSignature(modes[i], modes[j])
            if p: print(str(i) + ", " + str(j))
            if p: print(str(sig))
            for k in range(len(testBasis)):
                if polar(innerProduct(sig, testBasis[k]))[0] > 0.0001:
                    statesPresent.append(testBasis[k]) # (c, p)
            if p: print("### " + str(len(statesPresent)) + " ###")
            for state in statesPresent:
                if p: print(state)
            if p: print()
            largest = max(len(statesPresent), largest)

    spans = [False]*len(testBasis) 

    for k in range(len(testBasis)):
        for i in range(len(modes)):
            for j in range(i, len(modes)):
                sig = getSignature(modes[i], modes[j])
                if polar(innerProduct(sig, testBasis[k]))[0] > 0.0001:
                     spans[k] = True
    if p: print()
    if p: print(spans)                 
    if p: print("every Bell state present: " + str(sum(spans) == len(testBasis)))
    if p: print("largest: " + str(largest))
    if sum(spans) == len(testBasis):
        return(largest)
    else:
        return(999)

def multiTest():
    minSoFar = 1000
    subsets = list(combinations(bellBasis, 5))
    for i in range(len(subsets)):
        if i % 100 == 0:
            print(i)
        minSoFar = min(minSoFar, test(subsets[i]))
    return minSoFar






