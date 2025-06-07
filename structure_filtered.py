import argparse
import subprocess
from Bio import SeqIO

# Параметры по умолчанию
PROBE_LENGTH = 120
DEFAULT_DG_THRESHOLD = -9.0

def calculate_dg(sequence: str) -> float | None:
    """Вычисляет ΔG вторичной структуры с помощью RNAfold."""
    try:
        result = subprocess.run(
            ["RNAfold", "--noPS"],
            input=sequence,
            capture_output=True,
            text=True,
            check=True
        )
        output_lines = result.stdout.strip().split("\n")
        if len(output_lines) >= 2 and "(" in output_lines[1]:
            dg_str = output_lines[1].split()[-1].strip("()")
            return float(dg_str)
    except Exception as e:
        print(f"[!] Ошибка RNAfold: {e}")
    return None

def parse_args():
    parser = argparse.ArgumentParser(description="Фильтрация зондов по ΔG вторичной структуры (RNAfold)")
    parser.add_argument("input_fasta", help="Входной FASTA-файл")
    parser.add_argument("output_fasta", help="Выходной FASTA-файл")
    parser.add_argument("--dg-threshold", type=float, default=DEFAULT_DG_THRESHOLD,
                        help="Порог ΔG (по умолчанию: -9.0 ккал/моль)")
    return parser.parse_args()

def main():
    args = parse_args()

    input_file = args.input_fasta
    output_file = args.output_fasta
    dg_threshold = args.dg_threshold

    filtered_records = []
    total_checked = 0

    print(f"[•] Обработка {input_file}, порог ΔG: {dg_threshold} ккал/моль...")

    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq).upper()

        if len(seq) != PROBE_LENGTH:
            continue

        total_checked += 1
        dg = calculate_dg(seq)
        if dg is not None and dg > dg_threshold:
            record.description += f" ΔG={dg:.2f}"
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_file, "fasta")

    print(f"[✓] Сохранено {len(filtered_records)} зондов из {total_checked} (ΔG > {dg_threshold}) → {output_file}")

if __name__ == "__main__":
    main()
