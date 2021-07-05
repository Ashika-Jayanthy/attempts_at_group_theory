from Bio.Align.Applications import MuscleCommandline
import glob

indir = "./Data/sequence_by_family"
outdir = "./Data/alignments"

files = glob.glob(f"{indir}/*")
for family_fasta in files:
    name = family_fasta.split("/")[-1].split(".")[0]
    muscle_cline = MuscleCommandline(input=f"indir/{family_fasta}.fasta",out=f"outdir/{name}.aligned.fasta",)
    muscle_cline()
