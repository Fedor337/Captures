import argparse
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Фильтрация зондов по температуре плавления (Tm) из FASTA-файла."
    )
    parser.add_argument("input_fasta", help="Входной FASTA-файл с зондами")
    parser.add_argument("output_fasta", help="Выходной FASTA-файл для отфильтрованных зондов")
    parser.add_argument("--tm_min", type=float, default=65.0, help="Минимальная температура плавления (по умолчанию 65.0)")
    parser.add_argument("--tm_max", type=float, default=72.0, help="Максимальная температура плавления (по умолчанию 72.0)")
    return parser.parse_args()

def main():
    args = parse_args()
    filtered_probes = []

    try:
        records = list(SeqIO.parse(args.input_fasta, "fasta"))
    except Exception as e:
        sys.exit(f"[Ошибка] Не удалось прочитать FASTA: {e}")

    for record in records:
        sequence = str(record.seq).upper()
        if len(sequence) == 120 and set(sequence).issubset({"A", "T", "G", "C"}):
            tm = mt.Tm_NN(sequence)
            if args.tm_min <= tm <= args.tm_max:
                filtered_probes.append(record)

    SeqIO.write(filtered_probes, args.output_fasta, "fasta")
    print(f"[✓] {len(filtered_probes)} зондов записано в {args.output_fasta} (Tm от {args.tm_min} до {args.tm_max} °C)")

