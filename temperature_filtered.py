from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt

# Входной и выходной файлы
gc_filtered_input = "probes_gc_filtered.fa"
tm_filtered_output = "probes_gc__tm_filtered.fa"

# Пороговые значения температуры плавления
tm_min = 65.0
tm_max = 72.0

filtered_probes = []

for record in SeqIO.parse(gc_filtered_input, "fasta"):
    sequence = str(record.seq)
    
    # Фильтруем по длине
    if len(sequence) == 120:
        tm = mt.Tm_NN(sequence)  # Использует nearest-neighbor метод
        
        if tm_min <= tm <= tm_max:
            filtered_probes.append(record)

# Сохраняем отфильтрованные последовательности
SeqIO.write(filtered_probes, tm_filtered_output, "fasta")

print(f"{len(filtered_probes)} probes written to {tm_filtered_output} (length = 120 nt, Tm 65–72°C)")

