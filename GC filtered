from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq) * 100

input_file = "brca_probes.fa"
output_file = "probes_gc_filtered.fa"

filtered_probes = []

for record in SeqIO.parse(input_file, "fasta"):
    gc = gc_content(str(record.seq))
    if 40.0 <= gc <= 60.0:
        filtered_probes.append(record)

SeqIO.write(filtered_probes, output_file, "fasta")

print(f"{len(filtered_probes)} probes written to {output_file} (GC 40–60%)")
