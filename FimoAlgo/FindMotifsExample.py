#!/usr/bin/env python
# coding: utf-8

# <h3>Demonstrate the use of MotifsFinder</h3>

import numpy as np
import pandas as pd
from Bio import Align
from Motifstools.FindMotifs import FindMotifs


# Input examples. The files with the sequences must be in the Input directory in the format gene_name.fasta.
# The genes in sample_list_at/vv are the one that sould be homologues from literature.
# The task is one of the Phd thesis.

sample_input_list = ["meme_ara","meme_ara_1","meme_ara_2","meme_ara_3","meme_ara_4","meme_ara_5","meme_ara_6","meme_ara_7","meme_ara","meme_ara_1","meme_ara_2","meme_ara_3","meme_ara_4","meme_ara_5","meme_ara_6","meme_ara_7"]
sample_list_at = ["AT2G46680","AT4G26080","AT3G19290","AT1G08810","AT2G47160"]
sample_list_vv = ["VIT_15s0048g02870","VIT_11s0016g03180","VIT_18s0001g10450","VIT_08s0056g00800","VIT_17s0000g08530"]
sample_task_at = "AT4G28110"
sample_task_vv = ["VIT_12s0134g00570","VIT_19s0014g03820","VIT_00s0203g00070"]


# Used for alignment

at_cdna = pd.read_csv("Input/arabidopsis_cDNA.csv",sep="\t")
vv_cdna = pd.read_csv("Input/vitis_cDNA.csv",sep="\t")


# Used for some exmaples below

def align(seq1,seq2):
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -2
    aligner.match = 2
    aligner.mismatch = -3
    score = aligner.score(seq1, seq2)
    return score 


def jaccard_similarity(x,y):
    inter = len(set.intersection(set(x), set(y)))
    return inter / float(len(set(x)) + len(set(y)) - inter)


# Compute the motifs for two example lists (should be the expansion lists)

at = FindMotifs(sample_list_at)
vv = FindMotifs(sample_list_vv)


# Compute some values as testing.<br>
# sim_j: matrix 5x5 with only jaccard distances<br>
# sim_a: matrix 5x5 with only alignment<br>
# sim: matrix 5x5 with jaccard*alignmet<br>

sim_j = np.zeros((5,5))
for i in range(len(at)):
    for j in range(len(vv)):
        sim_j[i,j] = jaccard_similarity(at[i],vv[j])


sim_a = np.zeros((5,5))
for i in range(5):
    for j in range(5):
        a = at_cdna.loc[(at_cdna["a_thaliana_gene"] == sample_list_at[i]),"a_sequence"].iloc[0]
        v = vv_cdna.loc[(vv_cdna["vitis_gene"] == sample_list_vv[j]),"v_sequence"].iloc[0]
        sim_a[i,j] = align(a,v)


sim = np.zeros((5,5))
for i in range(len(at)):
    for j in range(len(vv)):
        a = at_cdna.loc[(at_cdna["a_thaliana_gene"] == sample_list_at[i]),"a_sequence"].iloc[0]
        v = vv_cdna.loc[(vv_cdna["vitis_gene"] == sample_list_vv[j]),"v_sequence"].iloc[0]
        sim[i,j] = jaccard_similarity(at[i],vv[j]) * align(a,v)


print("Jaccard")
print(sim_j)
print("Align")
print(sim_a)
print("Both")
print(sim)


# Compute the ratio between the max and the other of a line of the matrix to see if things get better from sim_j/a to sim

print("Jaccard")
print(sim_j[4,:] / max(sim_j[4,:]))
print("Align")
print(sim_a[4,:] / max(sim_a[4,:]))
print("Both")
print(sim[4,:] / max(sim[4,:]))


# Test for task

task_at = FindMotifs([sample_task_at])
task_vv = FindMotifs(sample_task_vv)


print(jaccard_similarity(task_at[0],task_vv[0]))
print(jaccard_similarity(task_at[0],task_vv[1]))
print(jaccard_similarity(task_at[0],task_vv[2]))

