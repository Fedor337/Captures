from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

# Параметры
gc_filtered_input = "probes_gc_filtered.fa"
tm_filtered_output = "probes_gc__tm_filtered.fa"

tm_min = 65.0
tm_max = 72.0

filtered_probes = []

for record in SeqIO.parse(gc_filtered_input, "fasta"):
    seq = str(record.seq)
    tm = mt.Tm_NN(seq)
    if tm_min <= tm <= tm_max:
        filtered_probes.append(record)

SeqIO.write(filtered_probes, tm_filtered_output, "fasta")

print(f"{len(filtered_probes)} probes written to {tm_filtered_output} (Tm 65–72°C)")
