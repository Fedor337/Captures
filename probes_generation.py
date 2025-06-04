from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import math

def make_probes(exon_seq, exon_id, probe_length=120, max_step=60):
    probes = []
    exon_len = len(exon_seq)

    if exon_len <= probe_length:
        # Если экзон короче зонда — создаём один зонд
        probes.append(SeqRecord(Seq(exon_seq), id=f"{exon_id}_probe_1_{exon_len}", description=""))
        return probes

    # Вычисляем необходимое количество зондов для покрытия всей последовательности
    num_probes = math.ceil((exon_len - probe_length) / max_step) + 1

    # Пересчитываем реальный шаг (не больше max_step)
    step = (exon_len - probe_length) / (num_probes - 1)
    step = int(math.floor(step))  # округляем вниз, чтобы не превысить max_step

    for i in range(num_probes):
        start = i * step
        end = start + probe_length
        if end > exon_len:
            # Смещаем последний зонд так, чтобы он полностью влез
            start = exon_len - probe_length
            end = exon_len
        probe_seq = exon_seq[start:end]
        probe_id = f"{exon_id}_probe_{start+1}_{end}"
        probes.append(SeqRecord(Seq(probe_seq), id=probe_id, description=""))
        if end == exon_len:
            break  # достигли конца — выходим
    return probes

# Основной блок
input_file = "brca_exons.fa"
output_file = "brca_probes.fa"
all_probes = []

for record in SeqIO.parse(input_file, "fasta"):
    exon_id = record.id
    exon_seq = str(record.seq)
    probes = make_probes(exon_seq, exon_id)
    all_probes.extend(probes)

# Сохраняем в FASTA
SeqIO.write(all_probes, output_file, "fasta")
print(f"{len(all_probes)} probes written to {output_file}")
