from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.lc_pattern import low_complexity_regions
import re

# Входной и выходной файлы
input_file = "brca_exons.fa"
output_file = "exons_filtered.fa"

# Параметры
PROBE_LENGTH = 120
HOMOPOLYMER_THRESHOLD = 5
TANDEM_MIN_REPEATS = 3
PALINDROME_MIN_LENGTH = 6

# --- Функции фильтрации ---
def has_homopolymers(seq: str) -> bool:
    """Проверка на гомополимеры (например, AAAAA)."""
    return bool(re.search(r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD
    ), seq))

def has_tandem_repeats(seq: str) -> bool:
    """Проверка на тандемные повторы (повтор мотивов от 2 до 5 нт)."""
    for size in range(2, 6):
        pattern = re.compile(r"((\w{%d}))\2{%d,}" % (size, TANDEM_MIN_REPEATS - 1))
        if pattern.search(seq):
            return True
    return False

def has_palindromes(seq_obj: Seq) -> bool:
    """Проверка на палиндромы длиной ≥ PALINDROME_MIN_LENGTH."""
    seq = str(seq_obj)
    for i in range(len(seq) - PALINDROME_MIN_LENGTH + 1):
        for l in range(PALINDROME_MIN_LENGTH, len(seq) - i + 1):
            fragment = seq[i:i + l]
            if fragment == str(Seq(fragment).reverse_complement()):
                return True
    return False

def has_low_complexity(seq_obj: Seq) -> bool:
    """Проверка на низкосложные регионы (DUST-подобный алгоритм)."""
    return len(low_complexity_regions(seq_obj)) > 0

# --- Основной процесс фильтрации ---
filtered_probes = []

for record in SeqIO.parse(input_file, "fasta"):
    seq = record.seq.upper()

    if len(seq) != PROBE_LENGTH:
        continue

    if has_homopolymers(str(seq)):
        continue
    if has_tandem_repeats(str(seq)):
        continue
    if has_palindromes(seq):
        continue
    if has_low_complexity(seq):
        continue

    filtered_probes.append(record)

# --- Сохранение результатов ---
SeqIO.write(filtered_probes, output_file, "fasta")

# --- Статистика ---
total = sum(1 for r in SeqIO.parse(input_file, "fasta") if len(r.seq) == PROBE_LENGTH)
print(f"Сохранено {len(filtered_probes)} из {total} зондов длиной 120 нт.")
