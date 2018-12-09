#!/usr/bin/env python
# coding: utf-8

import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from multiprocessing import Pool


motifs = "motifs_at.meme"
input_dir = "Input"


# Treshold for p value. Harcoded for now :)

p_treshold = 1e-3


# Given a list, eliminated duplicated by appending an incremental value

def resolve_duplicates(lst):
    d = {k:0 for k in list(set(lst))}
    for i in range(len(lst)):
        d[lst.iloc[i]]+=1
        lst.iloc[i] = lst.iloc[i] + "_r" + str(d[lst.iloc[i]])


# Given a gene name, find the motifs in the corresponding upstream sequence. 
# The file with the sequence must be inside the input_dir directory
# N.B. resolve_duplicated is commented out since it may yeld better results (?)

def compute_list(gene):
    try:
        res = subprocess.run(["fimo","--skip-matched-sequence", "--verbosity", "1", input_dir + "/" + motifs, input_dir + "/" + gene + ".fasta"],capture_output=True,check=True,text=True)
        result = pd.read_csv(StringIO(res.stdout), sep="\t", usecols=["motif_id","p-value"])
        result_cut = result.loc[result["p-value"] < p_treshold,"motif_id"]
        #resolve_duplicates(result_cut)
        return result_cut
    except subprocess.CalledProcessError as exp:
        return None


# Find the motifs present in each of the genes in gene_list.
# Return a Pandas Series with name_of_the_gene : list_of_motifs for each gene.

def FindMotifs(gene_list):
    with Pool(None) as p:
        results = p.map(compute_list, gene_list)
    return pd.Series(results, index=gene_list)

