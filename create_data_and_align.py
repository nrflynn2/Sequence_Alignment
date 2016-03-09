import numpy as np
from numpy import array
import sys
import random
from Bio import SeqIO
import subprocess

L = int(float(sys.argv[1]))

def generateSequence(L):
	seq = ""
	for i in range(L):
		seq += random.choice("ACTG")
	return seq

def mutateSequence(seq, L):
    mutations = int(L/10)
    newString = seq
    for mut in range(0, mutations):
        temp = ''
        position = random.randrange(0,len(newString))
        probability = int(2*random.random())
        if probability == 0:
            mutation = random.choice("ACTG")
            while mutation == newString[position]:
                mutation = random.choice("ACTG")
            temp = newString[:position] + mutation + newString[position+1:]
        elif probability == 1:
            temp = newString[:position] + newString[position+1:]
        newString = temp
    return newString

subject = generateSequence(L)
query = mutateSequence(subject, L)

#this function takes a name and sequence and returns a FASTA-formatted string
def format_fasta(name, sequence):
    fasta_string = '>' + name + "\n" + sequence + "\n"
    return fasta_string

#test the function by using it in a print statement
#print(format_fasta('sub', subject))
#print(format_fasta('que', query))

output_sub = open('sub.fasta', 'w')
output_que = open('que.fasta', 'w')
output_sub.write(format_fasta('sub', subject))
output_que.write(format_fasta('que', query))

#subprocess.call(["python", "align.py", "sub.fasta", "que.fasta", "subs.txt", "-500"])
subprocess.Popen("python align.py sub.fasta que.fasta subs.txt -500", shell=True)
