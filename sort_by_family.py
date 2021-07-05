from Bio import SeqIO

indir = "./Data"
outdir = "./Data/sequence_by_family"

families = []
for record in SeqIO.parse(f"{indir}/human_viruses.fasta", "fasta"):
    name = record.description.split("|")[1]
    if name not in families and len(name)>1:
        families.append(name.strip())

for family in families:
    records = []
    for record in SeqIO.parse(f"{indir}/human_viruses.fasta", "fasta"):
        name = record.description.split("|")[1]
        if name.strip() == family:
            records.append(record)
    SeqIO.write(records, f"{outdir}/{family}.fasta", "fasta")
