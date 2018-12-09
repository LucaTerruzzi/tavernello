#******************************************************************************
#                  FOF (Functional Orthologous Finder)                        #
#******************************************************************************
#
# This should be the main program to call 
from Bio import Align
import random as rnd #just to test time
import os

# example of SW call 
def align(seq1,seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -2
    aligner.match = 2
    aligner.mismatch = -3
    score = aligner.score(seq1, seq2)
    return score 

seq1 = ""
seq2 = ""
for i in range(1000): #tune it to test time dependency to length
    seq1 += "".join(rnd.sample("ACTG", 1))
    seq2 += "".join(rnd.sample("ACTG", 1))

start_time = os.times()[0]
print("Score: {}".format(align(seq1,seq2)))
print("Elapsed time (in seconds): {}".format(round(os.times()[0]-start_time, 3)))
