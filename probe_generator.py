from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import List, Union
from os import PathLike
import math


class ProbeGenerator:
    """
    Generates overlapping oligonucleotide probes from exon sequences.
    """

    def __init__(
        self,
        input_fasta: Union[str, PathLike] = "data/brca_exons.fa",
        output_fasta: Union[str, PathLike] = "data/brca_probes.fa",
        probe_length: int = 120,
        max_step: int = 60,
    ):
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.probe_length = probe_length
        self.max_step = max_step

    def make_probes(self, exon_seq: str, exon_id: str) -> List[SeqRecord]:
        """
        Generate overlapping probes from a single exon sequence.
        """
        probes = []
        exon_len = len(exon_seq)

        if exon_len <= self.probe_length:
            probes.append(SeqRecord(Seq(exon_seq), id=f"{exon_id}_probe_1_{exon_len}", description=""))
            return probes

        num_probes = math.ceil((exon_len - self.probe_length) / self.max_step) + 1
        step_size = int((exon_len - self.probe_length) / (num_probes - 1))

        for i in range(num_probes):
            start = i * step_size
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
        """
        Generate probes for all exons and write to output FASTA.
        """
        all_probes = []

        for record in SeqIO.parse(self.input_fasta, "fasta"):
            exon_id = record.id
            exon_seq = str(record.seq)
            exon_probes = self.make_probes(exon_seq, exon_id)
            all_probes.extend(exon_probes)

        SeqIO.write(all_probes, self.output_fasta, "fasta")
        print(f"[OK] {len(all_probes)} probes written to {self.output_fasta}")
