import sys
import random
from Bio import SeqIO
import numpy as np

fileOne = sys.argv[1]
fileTwo = sys.argv[2]
scoreMatrix = sys.argv[3]
gapPen = int(float(sys.argv[4]))

#seqOne = 'gattagc'
#seqTwo = 'gcattgc'

my_fileOne = open(fileOne)
seqOne = (my_fileOne.read()).rstrip("\n")
my_fileTwo = open(fileTwo)
seqTwo = (my_fileTwo.read()).rstrip("\n")

score = ((open(scoreMatrix)).read()).rstrip("\n")
sc = score.split('\t')
print(sc)

#Assigning values for match and mismatch as derived from the scoring matrix
AA = int(float(sc[5]))#Covers AA and TT
GG = int(float(sc[10])) #Covers GG and CC
AT = int(float(sc[17])) #Covers AT mismatch
AC = int(float(sc[6])) #Covers AC and GT mismatch
AG = int(float(sc[7])) #Covers AG and TC mismatch
GC = int(float(sc[11])) #Covers GC mismatch

# TO DO:
# Print out actual sequence alignment
# Implement input scoring matrix into scoring function

#=============================================================
# Alignment Parameters
#=============================================================
def generateSequence(L):
	seq = ""
	for i in range(L):
		seq += random.choice("ACTG")
	return seq

class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    #def __init__(self, match, mismatch, gap, gap_start=-2):
    def __init__(self, AA, GG, AT, AC, AG, GC, gap, gap_start=0):
        self.gap_start = gap_start
        self.gap = gap
        #self.match = match
        #self.mismatch = mismatch
        self.AA = AA
        self.GG = GG
        self.AT = AT
        self.AC = AC
        self.AG = AG
        self.GC = GC

    def matchchar(self, a,b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a=='a':
            if b=='a':
                return self.AA
            elif b=='t':
                return self.AT
            elif b=='c':
                return self.AC
            elif b=='g':
                return self.AG

        elif a=='t':
            if b=='t':
                return self.AA
            elif b=='a':
                return self.AT
            elif b=='c':
                return self.AG
            elif b=='g':
                return self.AC

        elif a=='g':
            if b=='g':
                return self.GG
            elif b=='t':
                return self.AC
            elif b=='c':
                return self.GC
            elif b=='a':
                return self.AG

        elif a=='c':
            if b=='c':
                return self.GG
            elif b=='a':
                return self.AC
            elif b=='g':
                return self.GC
            elif b=='t':
                return self.AG
        '''if a==b:
            return self.match
        else:
            return self.mismatch'''

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (
                self.match, self.mismatch, self.gap_start, self.gap
        )

#=============================================================
# Sequence Alignment
#=============================================================
def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]

def print_matrix(x, y, A):
    """Print the matrix with the (0,0) entry in the top left
    corner. Will label the rows by the sequence and add in the
    0-row if appropriate."""

    # decide whether there is a 0th row/column
    if len(x) == len(A):
        print "%5s" % (" "),
    else:
        print "%5s %5s" % (" ","*"),
        y = "*" + y

    # print the top row
    for c in x:
        print "%5s" % (c),
    print

    for j in xrange(len(A[0])):
        print "%5s" % (y[j]),
        for i in xrange(len(A)):
            print "%5.0f" % (A[i][j]),
        print

def local_align(x, y, score=ScoreParam(AA, GG, AT, AC, AG, GC, gapPen)):
    """Do a local alignment between x and y with the given scoring parameters.
    We assume we are MAXIMIZING."""

    # create a zero-filled matrix
    A = make_matrix(len(x) + 1, len(y) + 1)

    best = 0
    optloc = (0,0)

    alignOne = []
    alignTwo = []

    #trace = []

    # fill in A in the right order
    for i in xrange(1, len(x)+1):
        #iterTrace = []
        for j in xrange(1, len(y)+1):
            # the local alignment recurrance rule:
            A[i][j] = max(
               A[i][j-1] + score.gap,
               A[i-1][j] + score.gap,
               A[i-1][j-1] + score.matchchar(x[i-1], y[j-1]),
               #0
            )

            #Find traceback
            '''if A[i][j] ==  A[i][j-1] + score.gap:
                iterTrace.append((i, j-1))
            elif A[i][j] ==  A[i-1][j] + score.gap:
                iterTrace.append((i-1, j))
            elif A[i][j] ==  A[i-1][j-1] + score.matchchar(x[i-1], y[j-1]):
                iterTrace.append((i-1, j-1))'''

            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)
                alignOne.append(x[i-1])
                alignTwo.append(y[j-1])

        #trace.append(iterTrace)
    alOne = ''.join(alignOne)
    alTwo = ''.join(alignTwo)

    print "Scoring:", str(score)
    print "A matrix ="
    #print_matrix(x, y, A)
    print "Optimal Score =", best
    print "Max location in matrix =", optloc
    #for i in range(0, len(trace)):
        #print(trace[i])
    # return the opt score and the best location
    return best, optloc, alOne, alTwo

#result = local_align(seqOne, seqTwo, ScoreParam(gap=-2, match=2, mismatch=-1))
result = local_align(seqOne, seqTwo, ScoreParam(AA=AA, GG=GG, AT=AT, AC=AC, AG=AG, GC=GC, gap=gapPen, ))
print(result)
