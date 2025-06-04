import subprocess
from Bio import SeqIO

# Входной и выходной файлы
input_file = "probes_final_filtered.fa"
output_file = "probes_structure_filtered.fa"

# Порог свободной энергии ΔG (в ккал/моль)
DG_THRESHOLD = -9.0

def calculate_dg(sequence):
    """
    Вычисляет ΔG структуры с помощью RNAfold.
    Возвращает ΔG в ккал/моль или None при ошибке.
    """
    try:
        result = subprocess.run(
            ["RNAfold", "--noPS"],
            input=sequence,
            capture_output=True,
            text=True
        )
        output_lines = result.stdout.strip().split("\n")
        if len(output_lines) >= 2 and "(" in output_lines[1]:
            dg_str = output_lines[1].split()[-1].strip("()")
            return float(dg_str)
    except Exception as e:
        print(f"Ошибка при обработке последовательности: {e}")
    return None

# Чтение и фильтрация
filtered_records = []

print("Фильтрация зондов по вторичной структуре...")

for record in SeqIO.parse(input_file, "fasta"):
    seq = str(record.seq)
    dg = calculate_dg(seq)
    if dg is not None:
        if dg > DG_THRESHOLD:
            # Добавим ΔG в описание (опционально)
            record.description += f" ΔG={dg:.2f}"
            filtered_records.append(record)

# Сохранение отфильтрованных зондов
SeqIO.write(filtered_records, output_file, "fasta")

print(f"Готово! {len(filtered_records)} зондов сохранено в {output_file} (ΔG > {DG_THRESHOLD})")
