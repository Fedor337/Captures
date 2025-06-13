import argparse
from Bio import SeqIO

def gc_content(seq: str) -> float:
    """Вычисляет GC-содержание в процентах."""
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def parse_args():
    parser = argparse.ArgumentParser(
        description="Фильтрация зондов по GC-содержанию из FASTA-файла"
    )
    parser.add_argument("input_fasta", help="Входной FASTA-файл с зондами")
    parser.add_argument("output_fasta", help="Выходной FASTA-файл для отфильтрованных зондов")
    parser.add_argument("--gc_min", type=float, default=40.0, help="Минимальное GC-содержание (по умолчанию 40.0)")
    parser.add_argument("--gc_max", type=float, default=60.0, help="Максимальное GC-содержание (по умолчанию 60.0)")
    return parser.parse_args()

def main():
    args = parse_args()
    filtered_probes = []

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        sequence = str(record.seq).upper()
        if len(sequence) == 120 and set(sequence).issubset({"A", "T", "G", "C"}):
            gc = gc_content(sequence)
            if args.gc_min <= gc <= args.gc_max:
                filtered_probes.append(record)

    SeqIO.write(filtered_probes, args.output_fasta, "fasta")
    print(f"[✓] {len(filtered_probes)} зондов записано в {args.output_fasta} (GC от {args.gc_min}% до {args.gc_max}%)")



