import subprocess
import glob

indir = "./Data/hmm_sequences"
outdir = "./Data/hmm_sequences_msa"

files = glob.glob(f"{indir}/*")
for family_fasta in files:
    name = family_fasta.split(".")[0].strip()
    subprocess.run(f"python3 muscle.py --sequence {family_fasta} --email jayanthyashika@gmail.com --format fasta --outfile {outdir}/{name}", shell=True)
