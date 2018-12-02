#******************************************************************************
#                  FOF (Functional Orthologous Finder)                        #
#******************************************************************************
#
# This should be the main program to call 
from Accessors.SmithWaterman import SmithWaterman 

# example of SW call 
seq1 = "ATCTATCTATGGGATGCTAGCTTAGTGG"
seq2 = "AGCCTTAGCGTGATCGGGAGGCTTAGTG"
scores = [2,-3,2,4] # Match, Mismatch, Gap, GapExt
par_file = "./Accessors/params_data.csv"
print(SmithWaterman(seq1, seq2, scores, par_file))
