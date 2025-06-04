from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.lc_pattern import low_complexity_regions
import re

# Входной и выходной файлы
input_file = "brca_exons.fa"
output_file = "exons_filtered.fa"

# Параметры фильтрации
HOMOPOLYMER_THRESHOLD = 5
TANDEM_MIN_REPEATS = 3
PALINDROME_MIN_LENGTH = 6

def has_gomopolymers(seq):
    return bool(re.search(r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD
    ), seq))

def has_tandem_repeats(seq):
    for size in range(2, 6):  # Мотив 2–5 нт
        pattern = re.compile(r"((\w{%d}))\2{%d,}" % (size, TANDEM_MIN_REPEATS - 1))
        if pattern.search(seq):
            return True
    return False

def has_palindromes(seq_obj):
    seq = str(seq_obj)
    rev_comp = str(seq_obj.reverse_complement())
    for i in range(len(seq) - PALINDROME_MIN_LENGTH + 1):
        for l in range(PALINDROME_MIN_LENGTH, len(seq) - i + 1):
            fragment = seq[i:i + l]
            if fragment == Seq(fragment).reverse_complement():
                return True
    return False

def has_low_complexity(seq_obj):
    # Использует DUST-подобный алгоритм из Biopython
    regions = low_complexity_regions(seq_obj)
    return len(regions) > 0

# Основной цикл фильтрации
filtered_records = []
for record in SeqIO.parse(input_file, "fasta"):
    seq = record.seq.upper()

    if has_gomopolymers(str(seq)):
        continue
    if has_tandem_repeats(str(seq)):
        continue
    if has_palindromes(seq):
        continue
    if has_low_complexity(seq):
        continue

    filtered_records.append(record)

# Сохранение результатов
SeqIO.write(filtered_records, output_file, "fasta")
print(f"Сохранено {len(filtered_records)} из {sum(1 for _ in SeqIO.parse(input_file, 'fasta'))} экзонов.")
