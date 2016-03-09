import numpy as np
from numpy import array
import sys
import random
from Bio import SeqIO
import ast

A, C, G, T = 0, 1, 2, 3
int_to_char = {0:'A', 1:'C', 2:'G', 3:'T'}

scoreMatrix = sys.argv[3]
score = ((open(scoreMatrix)).read()).rstrip("\n")
sc = score.split('\t')

#Assigning values for match and mismatch as derived from the scoring matrix
AA = int(float(sc[5]))#Covers AA and TT
GG = int(float(sc[10])) #Covers GG and CC
AT = int(float(sc[17])) #Covers AT mismatch
AC = int(float(sc[6])) #Covers AC and GT mismatch
AG = int(float(sc[7])) #Covers AG and TC mismatch
GC = int(float(sc[11])) #Covers GC mismatch

best = 0
opt_loc = (0,0)

gap = int(float(sys.argv[4]))
#ACGT by ACGT
scoring = array([[AA,AC,AG,AT],
                 [AC,GG,GC,AG],
                 [AG,GC,GG,AC],
                 [AT,AG,AC,AA]])

#For referencing alignment program with homework problems
'''scoring = array([[2,-1,-1,-1],
                 [-1,2,-1,-1],
                 [-1,-1,2,-1],
                 [-1,-1,-1,2]])'''

fileOne = sys.argv[1]
fileTwo = sys.argv[2]
my_fileOne = open(fileOne)
seqOne = (my_fileOne.read()).rstrip("\n")
my_fileTwo = open(fileTwo)
seqTwo = (my_fileTwo.read()).rstrip("\n")

subject = (list(seqOne))[5:]
query = (list(seqTwo))[5:]
sub = []
que = []

def str_to_int(seq):
    result = []
    for char in seq:
        if char == 'A':
            result.append(0)
        elif char == 'C':
            result.append(1)
        elif char == 'G':
            result.append(2)
        elif char == 'T':
            result.append(3)
    return result

sub = str_to_int(subject)
que = str_to_int(query)

class AlignmentFinder(object):
    def __init__(self, seq1, seq2, best=0, opt_loc=(0,0)):
        self.seq1 = seq1
        self.seq2 = seq2
        self.D = None

    def find_gobal_alignment(self):
        self.D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self._compute_array()
        print self.D
        return self._traceback()

    def _compute_array(self):
        for i in xrange(self.seq1.size+1):
            self.D[i,0] = i*gap
        for j in xrange(self.seq2.size+1):
            self.D[0,j] = j*gap
        for i in xrange(1, self.seq1.size+1):
            for j in xrange(1, self.seq2.size+1):
                self.D[i,j] = max(  self.D[i-1, j-1] + self._get_score(i, j),
                                    self.D[i-1, j] + gap,
                                    self.D[i, j-1] + gap)

                # track the cell with the largest score
                if self.D[i,j] >= best:
                    self.best = self.D[i,j]
                    self.optloc = (i,j)

        print('The optimal alignment between given sequences has score ' + str(self.best) + '.')
        print('The matrix location of the optimal alignment score is ' + str(self.optloc) + '.')
    def _get_score(self, i, j):
        ''' To obtain the correct nucleotide in the sequence, we must
        substract 1 to the matrix index. '''
        return scoring[self.seq1[i-1], self.seq2[j-1]]

    def _get_aligned_pair(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '_'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '_'
        return (n1, n2)

    def _traceback(self):
        alignment= []
        i = self.seq1.size
        j = self.seq2.size
        while i >0 and j>0:
            if self.D[i-1, j-1] + self._get_score(i, j) == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, j))
                i -= 1
                j -= 1
            elif self.D[i-1, j] + gap == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, 0))
                i -= 1
            else:
                alignment.append(self._get_aligned_pair(0, j))
                j -= 1
        while i > 0:
            alignment.append(self._get_aligned_pair(i, 0))
            i -= 1
        while j > 0:
            alignment.append(self._get_aligned_pair(0, j))
            j -= 1
        alignment.reverse()
        return alignment

def print_sequences(pairs):
    top_seq = []
    bottom_seq = []
    for (b, t) in pairs:
        bottom_seq.append(b)
        top_seq.append(t)
    for n in top_seq:
        print n,
    print ' '
    for n in bottom_seq:
        print n,

if __name__ == "__main__":
    s1 = array(sub, dtype=np.int16)
    s2 = array(que, dtype=np.int16)
    aligner = AlignmentFinder(s1, s2)
    pairs = aligner.find_gobal_alignment()
    print_sequences(pairs)
