from Bio import SeqIO

def gc_content(seq: str) -> float:
    """Вычисляет GC-содержание в процентах."""
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

input_file = "brca_probes.fa"
output_file = "probes_gc_filtered.fa"

filtered_probes = []

for record in SeqIO.parse(input_file, "fasta"):
    sequence = str(record.seq)
    if len(sequence) == 120:
        gc = gc_content(sequence)
        if 40.0 <= gc <= 60.0:
            filtered_probes.append(record)

SeqIO.write(filtered_probes, output_file, "fasta")

print(f"{len(filtered_probes)} probes written to {output_file} (length = 120 nt, GC 40–60%)")
