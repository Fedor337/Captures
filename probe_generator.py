import argparse
from pathlib import Path
from typing import List, Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ProbeGenerator:
    def __init__(self, input_fasta: str, output_fasta: str, probe_length: int, step: int):
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.probe_length = probe_length
        self.step = step

        if not self.input_fasta.exists():
            raise FileNotFoundError(f"[Ошибка] Файл {self.input_fasta} не найден.")

    def make_probes(self, sequence: str, exon_id: str) -> List[SeqRecord]:
        probes = []
        for i in range(0, len(sequence) - self.probe_length + 1, self.step):
            probe_seq = sequence[i:i + self.probe_length]
            probe_id = f"{exon_id}_probe_{i+1}_{i + self.probe_length}"
            probes.append(SeqRecord(Seq(probe_seq), id=probe_id, description=""))
        return probes

    def generate_all(self):
        print(f"[🔍] Чтение FASTA: {self.input_fasta}")
        total_input = 0
        total_output = 0
        unique_sequences: Set[str] = set()
        all_probes: List[SeqRecord] = []

        for record in SeqIO.parse(self.input_fasta, "fasta"):
            total_input += 1
            exon_id = record.id
            sequence = str(record.seq).upper()

            if len(sequence) < self.probe_length:
                continue

            probes = self.make_probes(sequence, exon_id)
            for probe in probes:
                probe_seq_str = str(probe.seq)
                if probe_seq_str not in unique_sequences:
                    unique_sequences.add(probe_seq_str)
                    all_probes.append(probe)
                    total_output += 1

        if total_output == 0:
            print("[⚠] Ни одного уникального зонда не сгенерировано.")
        else:
            SeqIO.write(all_probes, self.output_fasta, "fasta")
            print(f"Сгенерировано {total_output} уникальных зондов.")
            print(f"Записано в файл: {self.output_fasta}")

        # Проверка существования выходного файла
        if not self.output_fasta.exists():
            print("Ошибка: файл зондов не создан.")
        else:
            print(f"Файл создан: {self.output_fasta.resolve()}")


def main():
    parser = argparse.ArgumentParser(description="Генерация уникальных зондов из FASTA-файла")
    parser.add_argument("input_fasta", help="Входной FASTA-файл с экзонами")
    parser.add_argument("output_fasta", help="Выходной FASTA-файл для зондов")
    parser.add_argument("--probe-length", type=int, default=120, help="Длина зонда (по умолчанию: 120)")
    parser.add_argument("--step", type=int, default=1, help="Шаг между зондами (по умолчанию: 1)")

    args = parser.parse_args()

    generator = ProbeGenerator(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        probe_length=args.probe_length,
        step=args.step
    )
    generator.generate_all()


if __name__ == "__main__":
    main()





