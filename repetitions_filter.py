import argparse
import re
import math
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

PROBE_LENGTH = 120

def has_homopolymers(seq: str, threshold: int) -> bool:
    """Проверка на гомополимеры длиной threshold и больше."""
    return bool(re.search(r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (
        threshold, threshold, threshold, threshold), seq))

def has_tandem_repeats(seq: str, min_repeats: int) -> bool:
    """Поиск тандемных повторов (2–5 нт)"""
    for size in range(2, 6):
        pattern = re.compile(r"((\w{%d}))\2{%d,}" % (size, min_repeats - 1))
        if pattern.search(seq):
            return True
    return False

def has_palindromes(seq_obj: Seq, min_length: int = 6) -> bool:
    """Проверка на палиндромные последовательности"""
    seq = str(seq_obj)
    for i in range(len(seq) - min_length + 1):
        for l in range(min_length, len(seq) - i + 1):
            fragment = seq[i:i + l]
            if fragment == str(Seq(fragment).reverse_complement()):
                return True
    return False

def shannon_entropy(seq: str) -> float:
    """Вычисляет энтропию Шеннона"""
    freq = Counter(seq)
    length = len(seq)
    return -sum((count / length) * math.log2(count / length) for count in freq.values())

def parse_args():
    parser = argparse.ArgumentParser(description="Фильтрация зондов от структурных артефактов")
    parser.add_argument("input_fasta", help="Входной FASTA-файл")
    parser.add_argument("output_fasta", help="Выходной FASTA-файл")
    parser.add_argument("--max-homopolymer", type=int, default=5, help="Максимальная допустимая длина гомополимера (по умолчанию 5)")
    parser.add_argument("--max-repeats", type=int, default=3, help="Минимальное число тандемных повторов для фильтрации (по умолчанию 3)")
    parser.add_argument("--min-entropy", type=float, default=1.8, help="Минимально допустимая энтропия Шеннона (по умолчанию 1.8)")
    parser.add_argument("--probe-length", type=int, default=120, help="Длина зонда (по умолчанию 120)")
    return parser.parse_args()

def main():
    args = parse_args()
    input_file = args.input_fasta
    output_file = args.output_fasta

    filtered_probes = []
    total = 0

    for record in SeqIO.parse(input_file, "fasta"):
        seq = record.seq.upper()
        if len(seq) != args.probe_length:
            continue
        total += 1

        if has_homopolymers(str(seq), args.max_homopolymer):
            continue
        if has_tandem_repeats(str(seq), args.max_repeats):
            continue
        if has_palindromes(seq):
            continue
        if shannon_entropy(str(seq)) < args.min_entropy:
            continue

        filtered_probes.append(record)

    SeqIO.write(filtered_probes, output_file, "fasta")
    print(f"[✓] Сохранено {len(filtered_probes)} из {total} зондов длиной {args.probe_length} нт")


