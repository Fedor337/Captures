import subprocess
from Bio import SeqIO

# Файлы
input_file = "Probes_final_filtered.fa"
output_file = "Probes_structure_filtered.fa"

# Порог по свободной энергии ΔG (в ккал/моль)
dg_threshold = -9.0
PROBE_LENGTH = 120

def calculate_dg(sequence):
    """
    Вычисляет ΔG вторичной структуры с помощью RNAfold.
    Возвращает ΔG (ккал/моль) или None в случае ошибки.
    """
    try:
        result = subprocess.run(
            ["RNAfold", "--noPS"],  # отключает генерацию файла с картинкой
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
        print(f"Ошибка при обработке последовательности: {e}")
    return None

# Основной процесс
filtered_records = []

print("Фильтрация зондов по вторичной структуре (ΔG)...")

for record in SeqIO.parse(input_file, "fasta"):
    seq = str(record.seq).upper()

    # Фильтрация по длине
    if len(seq) != PROBE_LENGTH:
        continue

    dg = calculate_dg(seq)
    if dg is not None and dg > dg_threshold:
        record.description += f" ΔG={dg:.2f}"
        filtered_records.append(record)

# Сохранение результатов
SeqIO.write(filtered_records, output_file, "fasta")

print(f"Готово! {len(filtered_records)} зондов сохранено в {output_file} (ΔG > {dg_threshold})")
