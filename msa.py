import subprocess
import glob

indir = "./Data/sequence_by_family"
outdir = "./Data/alignments"

files = glob.glob(f"{indir}/*")
for family_fasta in files:
    name = family_fasta.split("/")[-1].split(".")[0]
    subprocess.run(f"python3 muscle.py --sequence {family_fasta} --email jayanthyashika@gmail.com --format fasta --outfile {outdir}/{name}.aligned.fasta", shell=True)
