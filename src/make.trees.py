#!/usr/bin/env python3

# Prune all 1a sequences which are from the same patient (by a distance
# cutoff).

from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from io import StringIO
import itertools
import subprocess
import tempfile
import os
import sys
import csv
import random

def makedir(str):
    try:
        os.mkdir(str)
    except OSError:
        0

FASTTREE = "fasttree2"

align_file = sys.argv[1]
CUTOFF_DIR = sys.argv[2]

makedir(CUTOFF_DIR)

for rep in range(1, 1001):
    random.seed(rep)
    
    print(rep)
    
    out_file = os.path.join(CUTOFF_DIR, os.path.basename(align_file).replace(".fasta", "_" + str(rep) + ".nwk"))
    
    if os.path.exists(out_file):
        print("File exists")
        continue
        
    sys.stdout.flush()
    
    queries = list(SeqIO.parse(align_file, "fasta"))
    
    
    if (len(queries) < 3):
        continue
         
    swap = [i for i in range(len(queries[0]))]
    for i in range(len(queries[0])):
        r = random.randrange(i, len(queries[0]))
        tmp = swap[i]
        swap[i] = swap[r]
        swap[r] = tmp

    for i in range(len(queries)):
    	queries_old = queries[i]
    	for j in range(len(queries_old)):
    	    queries[i] = queries[i][:j] + queries_old[swap[j]] + queries[i][(j+1):]
    
    stdin = StringIO()
    SeqIO.write(queries, stdin, "fasta")

    fasttree = subprocess.Popen([FASTTREE, "-nt", "-gtr", "-nosupport", "-seed", str(rep)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, stdin=subprocess.PIPE, universal_newlines=True)
    tree_str = StringIO(fasttree.communicate(stdin.getvalue())[0])
    tree = Phylo.read(tree_str, "newick")
    Phylo.write(tree, out_file, "newick")
