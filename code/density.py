import numpy as np
import math


#     sample vectors written in std basis
#     00 01 02 03 10 11 12 13 20 21 22 23 30 31 32 33
v1 = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
v2 = [0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0]
v3 = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0]
v4 = [0, 0, 0, 1,-1j,0, 0, 0, 0,-1, 0, 0, 0, 0,1j, 0]

#      hyperentangled stuff written in bell basis (top bottom)
#      00 10 20 30 01 11 21 31 02 12 22 32 03 13 23 33
h1  = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
h2  = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
h3  = [0,1-1j,0,1+1j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
h4  = [0,1+1j,0,1-1j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
h5  = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,-1, 0]
h6  = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,-1, 0, 1, 0] 
h7  = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,-1j, 0,1j]
h8  = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,1j, 0,-1j]
h9  = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
h10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
h11 = [0, 0, 0, 0, 0, 0, 0, 0,0,1-1j,0,1+1j, 0, 0, 0, 0]
h12 = [0, 0, 0, 0, 0, 0, 0, 0,0,1+1j,0,1-1j, 0, 0, 0, 0]
h13 = [0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 1, 0, 1, 0] 
h14 = [0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0,-1, 0,-1, 0] 
h15 = [0, 0, 0, 0, 0,1j, 0,-1j,0, 0, 0, 0, 0,-1, 0,-1]
h16 = [0, 0, 0, 0, 0,1j, 0,-1j,0, 0, 0, 0, 0, 1, 0, 1]

## (the important 7) hyperentangled states in the std basis
#        00 01 02 03 10 11 12 13 20 21 22 23 30 31 32 33
hstd1 = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
hstd2 = [0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0]
hstd3 = [0, 0, 1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0,-1, 0, 0]
hstd4 = [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0]
hstd5 = [0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0]
hstd6 = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0]
hstd7 = [0, 0, 0, 1, 0, 0, 1, 0, 0,-1, 0, 0,-1, 0, 0, 0]

##################### test a hyperentangled -> bell transformation #####################

lynnsTransform =   [[1,1,0,0],   # 0 -> 0 + 1
                   [1,-1,0,0],   # 1 -> 0 - 1
                   [0,0,1,1],    # 2 -> 2 + 3
                   [0,0,1,-1]]   # 3 -> 2 - 3

tommyTransform1 =  [[1,1,1,1],      # 0 -> 0 +  1 +  2 +  3
                   [1,1j,-1,-1j],   # 1 -> 0 + i1 -  2 - i3
                   [1,-1,1,-1],     # 2 -> 0 -  1 +  2 -  3
                   [1,-1j,-1,1j]]   # 3 -> 0 - i1 -  2 + i3

tommyTransform2 =  [[1,1,1,1],      # 0 -> 0 + 1 + 2 + 3
                   [1,1,-1,-1],     # 1 -> 0 + 1 - 2 - 3
                   [1,-1,-1,1],     # 2 -> 0 - 1 - 2 + 3
                   [1,-1,1,-1]]     # 3 -> 0 - 1 + 2 - 3





def transform(hyp):

    oneKetTrans = tommyTransform2  # tells me how each single-particle ket transforms
    
    sumVec = [0] * 16

    # go through all the two-particles kets in hyp...
    for i in range(len(hyp)):
        if hyp[i] == 0: continue
        left = (int) (i/4)  # state of the left particle
        right = i % 4  # state of the right particle

        # transform the kets, and take their product
        prod = tensorProd(oneKetTrans[left],oneKetTrans[right])

        # add this to our final state
        for j in range(len(sumVec)):
            sumVec[j] += hyp[i] * prod[j]
    
    return sumVec

def tensorProd(left,right):
    ''' takes 2 single-particle states (BOTH vectors of len d) <-- pls don't fight me on this
        returns their product (vector of len d^2)
        ex: left:  |0> - |1> = [1,-1,0]
            right: |2>       = [0,0,1]
            out:   |02> - |12> = [0,0,1,0,0,-1,0,0,0]  '''
    d = len(left)
    prod = []
    
    for i in range(d):
        for j in range(d):
            prod.append(left[i]*right[j])
    
    return prod




##################### hyperentangled -> qu4it ############################

phi0 = [[1,0,0],[1,1,1]]
phi1 = [[1,0,0],[-1,1,1]]
psi0 = [[1,1,0],[1,0,1]]
psi1 = [[1,1,0],[-1,0,1]]

def bellCross(b1,b2):
    ''' converts a cross product of two qubit bells states to the corresponding superposition of qu4it bell states'''
    ''' b1, b2 are bell states, represented as [[coefficient, L, R],[coefficient, L, R]] '''
    prod = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    # foil it out
    for ket1 in b1:
        for ket2 in b2:
            newCoeff = ket1[0] * ket2[0]
            newL = 2*ket1[1] + ket2[1]
            newR = 2*ket1[2] + ket2[2]
            # so the ket we've arrived at, in the d=4 std basis, is |newL, newR>
            # now we just add this ket to the product
            pos = 4*newL + newR
            prod[pos] = prod[pos]+newCoeff
    return stdToBell(prod)

##################### basis conversion ############################

M = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,-1, 0, 1, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0,1+1j, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0,1-1j, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
     
     


# takes you from the bell basis to the std basis
# this is also the conjugate transpose of Pbellstd
# should be halved!
Pstdbell = [[1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,1j,-1,-1j],
            [1,1j,-1,-1j,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1,1j,-1,-1j,0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1,1j,-1,-1j,0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1],
            [1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1,-1j,-1,1j,0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1,-1j,-1,1j,0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,-1j,-1,1j],
            [1,-1j,-1,1j,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]


def bellToStd(b):
    ''' converts a bell basis vector to the standard basis '''
    return np.matmul(Pstdbell,b)

def stdToBell(s):
    ''' converts a std basis vector to the bell basis '''
    v = []
    for n in range(len(Pstdbell)):
        sum = 0
        for rowNum in range(len(Pstdbell)):
            sum += np.conjugate(Pstdbell[rowNum][n]) * s[rowNum]
        v += [sum]
    return v

##################### density matrix ############################

def trL(v):
    ''' gives the left trace of the density matrix corresponding to the state v
        v should be in the d=3 standard basis
        returns as a list containing the rows as lists'''
    ''' this is how you check if its fully entangled '''
    mat = []
    for r in range(4):
        row = []
        for c in range(4):
            el = 0
            for k in range(4):
                el += v[4*k+r] * np.conjugate(v[4*k+c])
            row += [el]
        mat += [row]
    return mat

################## prettying up outputs #########################

def prettyMatrix(m):
    ''' makes a string of the matrix m with complex entries, but pretty '''
    ''' this is useful for when ur using the trL density stuff '''
    out = ""
    for row in m:
        for z in row:
            add = ""
            if (z==0): add += '0'
            elif (np.imag(z)==0): add += str(np.real(z))
            elif (np.real(z)==0): add += str(np.imag(z))
            else: add += str(z)
            while len(add) < 5:
                add += " "
            out += add
        out += '\n'
    return out

def prettyd4(v):
    ''' makes a string of v, a vector of len 16 representing a d=4 2-particle state,
        in the std basis '''
    out = ""
    # make it extra pretty by dividing out common factors
    gcf = 0
    for i in range(len(v)):
        gcf = math.gcd(gcf, int(v[i].real))
        gcf = math.gcd(gcf, int(v[i].imag))
        #print(gcf)
    
    for i in range(len(v)):
        if v[i] == 0: continue
        right = i % 4
        left = (int) (i/4)
        coeffString = ""
        
        coeff = v[i]
        if gcf != 0: coeff = v[i]/gcf

        if coeff == 1: coeffString = " + "
        elif coeff == -1: coeffString = " - "
        else: coeffString = " + "+str(coeff) 

        out += coeffString+"|"+str(left)+str(right)+">"
    return out

##################### machine learning ############################

def doMachineLearningStep(transform, inputBasis, goalBasis):
    outputBasis = list(map(lambda x: np.matmul(transform, x), inputBasis))
    print(outputBasis)
    return(basisMeasure(inputBasis, goalBasis))

def basisMeasure(basis1, basis2):
    '''takes in two sets of vectors {v_1,..., v_n} and
                                    {u_1,..., u_n}
       returns the norm of the tuple of inner products, i.e.,
       sqrt(|<v_1, u_1>|^2 + ... + |<v_n, u_n>|^2)
       '''
    dim = len(basis1)

    # normalize the bases
    normBasis1 = list(map(lambda x: x / (np.linalg.norm(x)), basis1))
    normBasis2 = list(map(lambda x: x / (np.linalg.norm(x)), basis2))

    print(normBasis1)

    # add their inner products
    counter = []
    for i in range(dim):
        counter.append((np.inner(normBasis1[i], normBasis2[i]))**2)
    return(math.sqrt(sum(counter)/dim))


        