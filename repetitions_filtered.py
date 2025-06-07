import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.lc_pattern import low_complexity_regions

# Параметры
PROBE_LENGTH = 120
HOMOPOLYMER_THRESHOLD = 5
TANDEM_MIN_REPEATS = 3
PALINDROME_MIN_LENGTH = 6

def has_homopolymers(seq: str) -> bool:
    return bool(re.search(r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD,
        HOMOPOLYMER_THRESHOLD
    ), seq))

def has_tandem_repeats(seq: str) -> bool:
    for size in range(2, 6):
        pattern = re.compile(r"((\w{%d}))\2{%d,}" % (size, TANDEM_MIN_REPEATS - 1))
        if pattern.search(seq):
            return True
    return False

def has_palindromes(seq_obj: Seq) -> bool:
    seq = str(seq_obj)
    for i in range(len(seq) - PALINDROME_MIN_LENGTH + 1):
        for l in range(PALINDROME_MIN_LENGTH, len(seq) - i + 1):
            fragment = seq[i:i + l]
            if fragment == str(Seq(fragment).reverse_complement()):
                return True
    return False

def has_low_complexity(seq_obj: Seq) -> bool:
    return len(low_complexity_regions(seq_obj)) > 0

def parse_args():
    parser = argparse.ArgumentParser(description="Фильтрация зондов от структурных артефактов")
    parser.add_argument("input_fasta", help="Входной FASTA-файл")
    parser.add_argument("output_fasta", help="Выходной FASTA-файл")
    return parser.parse_args()

def main():
    args = parse_args()
    input_file = args.input_fasta
    output_file = args.output_fasta

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

    SeqIO.write(filtered_probes, output_file, "fasta")

    total = sum(1 for r in SeqIO.parse(input_file, "fasta") if len(r.seq) == PROBE_LENGTH)
    print(f"[✓] Сохранено {len(filtered_probes)} из {total} зондов длиной {PROBE_LENGTH} нт")

if __name__ == "__main__":
    main()
