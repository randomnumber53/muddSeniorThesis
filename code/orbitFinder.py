import copy
# start with initial bell state
# first array is phase
# second array is correlation
# an array [[a,b],[c,d]]
#           phase: ap + bc
#           correlation: cp + dc
# where p is the initial phase and c is the initial correlation

#################### HELPER METHODS ########################

# addState will add the bell state to nextStates
# if we haven't seen it before
def addState(state):
    if state not in totalStates:
        nextStates.append(state)
        totalStates.append(state)

# if the state is "new" (doesn't appear in totalStates)
# and add it to the new current states
# applySc will shift the correlation by phase
# and return the state
def applySc(state):
    state[1][0] += state[0][0]
    state[1][1] += state[0][1]
    # mod by d to prevent int overflow
    state[1][0] %= d
    state[1][1] %= d
    return state

# applySp will shift the phase by correlation
# and return the state
def applySp(state):
    state[0][0] += state[1][0]
    state[0][1] += state[1][1]
    # mod by d to prevent int overflow
    state[0][0] %= d
    state[0][1] %= d
    return state

#################################################

# We start with the initial state
currentStates = [[[1,0],[0,1]]]
d = 3

# states we have seen before
totalStates = [[[1,0],[0,1]]]

# states of the next generation
nextStates = []

# starting with the initial bell state,
# we will apply Sc (shift correlation by phase)
# and Sp (shift phase by correlation)
#   we stop if we reach:
#   a bell state we reached before mod d
# We will keep going until currentStates is empty

while len(currentStates) > 0:
    nextStates = []
    print(totalStates)
    for bellState in currentStates:
        # copy the bell state because python copies by reference :(
        copyState1 = copy.deepcopy(bellState)
        copyState2 = copy.deepcopy(bellState)
        # apply the group operations
        phaseShift = applySp(copyState1)
        corrShift = applySc(copyState2)
        # add the resulting states to totalStates if they are new
        addState(phaseShift)
        addState(corrShift)
    currentStates = nextStates

print(len(totalStates))


