'''
This script implements the Needleman-Wunsch algorithm as described on wikipedia:
 https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

'''

import numpy as np
gap = -6
endGap = 0

def read_sequences(filename):
    with open(filename, 'r') as f:
        s1,s2 = '', ''
        for i, line in enumerate(f):
            if i == 1 :
                s1 = line.rstrip()
            if i == 3:
                s2 = line.rstrip()
        return s1, s2

def read_pam(filename):
    pam = {}
    with open(filename, 'r') as f:
        firstLine = []
        for i, line in enumerate(f):
            line = line.rstrip()
            if i == 0 :
                firstLine = line.split('\t')
            else :
                tokens = line.split('\t')
                for j in range(len(tokens)):
                    if j > 0:
                        key = tokens[0] + firstLine[j]
                        pam[key] = tokens[j]
    return pam

def align(s1, s2):
    m = initializeMatrix(s1, s2)
    for i in range(1, m.shape[0]) :
        for j in range (1,m.shape[1]):
            m[i,j] = score(m,i,j)
    print m
    return m


def score(m, i, j):
    diagonalScore = m[i-1, j-1] + pamScoreAtIndex(i,j)
    leftScore = m[i, j-1] + gap
    topScore = m[i-1, j] + gap
    return max(diagonalScore, leftScore, topScore)

def pamScoreAtIndex(i,j):
    # we substract 1 from indexes because in the matrix we prefilled the first columns
    # and return 0 for the first element which doesn't have a predecessor
    if i == 0 or j == 0:
        return 0
    key = s1[i-1] + s2[j-1]
    return int(pam[key])

def initializeMatrix (s1, s2):
    m = np.zeros((len(s1)+1, len(s2)+1))
    for i in range(len(s1) + 1):
        m[i, 0 ] = i * endGap

    for i in range(len(s2) + 1):
        m[0,i] = i * endGap
    return m

def alignedSequences1(m) :
    r1, r2 = '', ''
    i, j =  maxScoreCoordinates(m)
    while i > 0 and j > 0 :
        diagonalScore = m[i-1, j-1]
        leftScore = m[i, j-1]
        topScore = m[i-1,j]
        pamScore = pamScoreAtIndex(i,j)
        currentScore = m[i, j]
        if currentScore == (diagonalScore + pamScore):
            r1 = s1[i-1] + r1
            r2 = s2[j-1] + r2

            i = decrementIndex(i)
            j = decrementIndex(j)
        if currentScore == leftScore + gap :
            r1 = '-' + r1
            j = decrementIndex(j)
        if currentScore == topScore + gap:
            r1 = s1[i-1] + r1
            r2 = '-' + r2
            i = decrementIndex(i)

    return r1, r2

def alignedSequences(m) :
    r1, r2 = '', ''

    i, j =  maxScoreCoordinates(m)

    while len(r1) < len(s1) :
        diagonalScore = m[i-1, j-1]
        leftScore = m[i, j-1]
        topScore = m[i-1,j]
        pamScore = pamScoreAtIndex(i,j)
        currentScore = m[i, j]

        if currentScore == (diagonalScore + pamScore):
            #diagonal
            r1 = s1[i-1] + r1
            r2 = s2[j-1] + r2
            i -= 1
            j -= 1

        if currentScore == leftScore + gap :
            #left
            r1 = '-' + r1
            j -= 1

        if currentScore == topScore + gap:
            #top
            r1 = s1[i-1] + r1
            r2 = '-' + r2
            i -=1

        #The following 2 cases occur when we reached the border of the matrix and on one sequence we have to fill in
        #with gaps and on the other one with the remaining letters obtaining values such as : HEAGAWGHEE ,---PAW-HEA
        if j == 0 :
            while i > j :
                r2 = '-' + r2
                r1 = s1[i-1] + r1
                i -=1

        if i == 0 :
            while j > i :
                r1 = '-' + r1
                r2 = s2[j-1] + r2
                j -=1

    return r1, r2

def decrementIndex(index):
    if index > 0:
        return index -1
    return index

def alignmentScore(m):
    return m.max()

def maxScoreCoordinates(m):
    maxScore = m.max()
    #loop in reverse order to find biggest values
    #Max score is not necessarily at m[n,n]
    for i in range(m.shape[0] -1 , 1, -1) :
        for j in range (m.shape[1]-1, 1,  -1):
            if m[i,j] == maxScore:
                return i, j

def alignmentScore(m):
    return m.max()




#s1 = 'HEAGAWGHEE'
#s2 = 'PAWHEAE'

s1,s2 = read_sequences('sequences.txt')
pam = read_pam('pam250.tab')
m = align(s1,s2)
print 'Score ' , alignmentScore(m)
print alignedSequences(m)
