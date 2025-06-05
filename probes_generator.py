from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
from typing import List


class ProbeGenerator:
    def __init__(self,
                 input_fasta: str | Path = "data/brca_exons.fa",
                 output_fasta: str | Path = "data/brca_probes.fa",
                 probe_length: int = 120,
                 max_step: int = 60):
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.probe_length = probe_length
        self.max_step = max_step

    def make_probes(self, exon_seq: str, exon_id: str) -> List[SeqRecord]:
        probes = []
        exon_len = len(exon_seq)

        if exon_len <= self.probe_length:
            probes.append(SeqRecord(Seq(exon_seq), id=f"{exon_id}_probe_1_{exon_len}", description=""))
            return probes

        num_probes = math.ceil((exon_len - self.probe_length) / self.max_step) + 1
        step = (exon_len - self.probe_length) / (num_probes - 1)
        step = int(math.floor(step))

        for i in range(num_probes):
            start = i * step
            end = start + self.probe_length
            if end > exon_len:
                start = exon_len - self.probe_length
                end = exon_len
            probe_seq = exon_seq[start:end]
            probe_id = f"{exon_id}_probe_{start+1}_{end}"
            probes.append(SeqRecord(Seq(probe_seq), id=probe_id, description=""))
            if end == exon_len:
                break
        return probes

    def generate_all(self) -> None:
        all_probes = []

        for record in SeqIO.parse(self.input_fasta, "fasta"):
            exon_id = record.id
            exon_seq = str(record.seq)
            exon_probes = self.make_probes(exon_seq, exon_id)
            all_probes.extend(exon_probes)

        SeqIO.write(all_probes, self.output_fasta, "fasta")
        print(f"[OK] {len(all_probes)} probes written to {self.output_fasta}")
