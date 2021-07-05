from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import SeqIO
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



n_seqs = 5
length = 5
def DNA(l):
    return ''.join(random.choice('CGTA') for _ in range(l))

records = []

for n in range(n_seqs):
    records.append(SeqRecord(Seq(DNA(length)),id=str(n)))

with open("test.fasta", "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")



muscle_cline = MuscleCommandline(input="test.fasta",out="outest.fasta",)
muscle_cline()
