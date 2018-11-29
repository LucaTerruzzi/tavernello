#*****************************#
# Smith Waterman Implentation #
#*****************************#

# Smith Waterman giving the Bit score, used Blast estimates for Lambda and K.
# Notice that bit score = log_2(q_ij/(p_i*p_j)) i.e lod score.
# Or, the bit-score is the required size of a sequence database 
# in which the current match could be found just by chance.
# It is preferred to just Alginement score since it's normalized and scaled 
# in case of different scores used. 
# With this implementation is possible to compute the p-value, check 
# https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html

import math 
            
# Accessory Functions ----------------------------------------------------------

def getParams(par_file, scores):
        "Get Lambda and K for Bit score computation from csv file"
        fh = open(par_file)
        found = False
        l = 0
        k = 0
        scores = [str(i) for i in scores]
        for line in fh:
            line = line.strip().split(",")
            if line[0:4] == scores:
                found = True 
                l = float(line[5])
                k = float(line[6])
        
        if not found:
            raise(Exception("Lamda and k Not available for such scores!"))

        return l,k

def findMax(matrix):
    """ Finds max in score matrix, if more than one,
        stores them in a list"""
    max_so_far = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > max_so_far:
                max_so_far = matrix[i][j]
                
    return max_so_far

# Matrix Construction ------------------------------------------------------------------------------

def checkEvents(pos, matrix):
    """ checks the max val for gap or mm from current pos.
        Stores data in a list composed by three vals:
        1) Score inside of matrix
        2) Position considered
        3) Type of event """

    g1 = [matrix[pos[0]-1][pos[1]], [pos[0]-1, pos[1]], "G"]
    g2 = [matrix[pos[0]][pos[1]-1], [pos[0], pos[1]-1], "G"]
    mm = [matrix[pos[0]-1][pos[1]-1], [pos[0]-1, pos[1]-1], "MM"]

    dict_moves = {g1[0]:g1[1:], g2[0]:g2[1:], mm[0]:mm[1:]}
    max_move = max(dict_moves)
    event = dict_moves[max_move][1]
    
    return event

def ConstructMatrix(str1, str2, M, MM, G, GE):
    """Constructing score matrix with M, MM, GG as penalties
    given in input. str1, str2 two sequence to align."""

    ls1 = len(str1)
    ls2 = len(str2)
    matrix = [[0]*(ls2+1) for i in range(ls1+1)]

    # filling the matrix
    for i in range(1,ls1+1): 
        for j in range(1,ls2+1):
            
            if str1[i-1] == str2[j-1]:
                matrix[i][j] = (matrix[i-1][j-1] + M)
                
            else:
                gaps = [[i-1,j],[i,j-1]]
                events = [checkEvents(g, matrix) for g in gaps]
                mm = matrix[i-1][j-1] + MM
                
                for event in events:
                    if event == "G":
                        g1 = matrix[i-1][j] - GE
                        g2 = matrix[i][j-1] - GE
                    else:
                        g1 = matrix[i-1][j] - G 
                        g2 = matrix[i][j-1] - G 
                
                matrix[i][j] = max(g1,g2,mm,0)

    return matrix
    
# Program Call -------------------------------------------------------------------------------------

def SmithWaterman(str1, str2, scores, par_file):
    """ Returns the alignemnt score in terms of bit score """
    match = scores[0]
    mismatch = scores[1]
    gap = scores[2]
    gap_ext = scores[3]
    matrix = ConstructMatrix(str1, str2, match, mismatch, gap, gap_ext) 
    score = findMax(matrix)
    l , k = getParams(par_file, scores)
    norm_score = (score*l - math.log(k))/math.log(2)
    
    return round(norm_score, 2)

