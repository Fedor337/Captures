import argparse
from pathlib import Path
from typing import List, Union, Set
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ProbeGenerator:
    """
    Generates overlapping probes (sliding window) from exon sequences.
    Keeps only unique probe sequences.
    """

    def __init__(
        self,
        input_fasta: Union[str, Path],
        output_fasta: Union[str, Path],
        probe_length: int = 120,
        step: int = 1,
    ):
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.probe_length = probe_length
        self.step = step

    def make_probes(self, exon_seq: str, exon_id: str) -> List[SeqRecord]:
        """Generate all overlapping probes from a single exon sequence."""
        probes = []
        exon_len = len(exon_seq)

        if exon_len < self.probe_length:
            return probes

        for start in range(0, exon_len - self.probe_length + 1, self.step):
            end = start + self.probe_length
            probe_seq = exon_seq[start:end]
            probe_id = f"{exon_id}_probe_{start+1}_{end}"
            probes.append(SeqRecord(Seq(probe_seq), id=probe_id, description=""))

        return probes

    def generate_all(self) -> None:
        """Generate unique probes for all exons and write to output FASTA."""
        unique_probes: List[SeqRecord] = []
        seen_seqs: Set[str] = set()

        for record in SeqIO.parse(self.input_fasta, "fasta"):
            exon_id = record.id
            exon_seq = str(record.seq).upper()
            exon_probes = self.make_probes(exon_seq, exon_id)

            for probe in exon_probes:
                seq_str = str(probe.seq)
                if seq_str not in seen_seqs:
                    seen_seqs.add(seq_str)
                    unique_probes.append(probe)

        SeqIO.write(unique_probes, self.output_fasta, "fasta")
        print(f"[✓] {len(unique_probes)} unique probes written to {self.output_fasta}")


def main():
    parser = argparse.ArgumentParser(description="Generate overlapping unique probes from exon sequences.")
    parser.add_argument("input_fasta", help="Input FASTA file with exon sequences")
    parser.add_argument("output_fasta", help="Output FASTA file to save unique probes")
    parser.add_argument("--probe-length", type=int, default=120, help="Length of each probe (default: 120)")
    parser.add_argument("--step", type=int, default=1, help="Step size between probes (default: 1 for 99% overlap)")

    args = parser.parse_args()

    pg = ProbeGenerator(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        probe_length=args.probe_length,
        step=args.step
    )
    pg.generate_all()


if __name__ == "__main__":
    main()




