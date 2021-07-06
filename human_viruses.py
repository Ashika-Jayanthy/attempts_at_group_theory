from Bio import pairwise2
import pyhmmer
from Bio import SeqIO
from core import *

def pairwise_alignment(s1,s2):
    score = pairwise2.align.globalxx("ACCGT", "ACG", score_only=True, one_alignment_only=True)
    return score
