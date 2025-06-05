from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import List
import re
import subprocess
import tempfile

class ProbeFilterPipeline:
    """
    Applies sequential filtering steps to a set of probes in FASTA format.

    Filters:
    - GC content (default: 40–60%)
    - Melting temperature (Tm, default: 65–72°C)
    - Repeats: homopolymers, di-/tri-nucleotide patterns
    - Secondary structure: RNAfold ΔG (optional)
    """
    def __init__(self,
                 input_fasta: str | Path,
                 output_fasta: str | Path,
                 gc_min: float = 40.0,
                 gc_max: float = 60.0,
                 tm_min: float = 65.0,
                 tm_max: float = 72.0,
                 allow_repeats: bool = False,
                 structure_filter: bool = False,
                 dg_threshold: float = -9.0):
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.allow_repeats = allow_repeats
        self.structure_filter = structure_filter
        self.dg_threshold = dg_threshold

    def filter_by_gc(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Retain probes with GC content within the specified range.
        """
        return [r for r in probes if self.gc_min <= GC(r.seq) <= self.gc_max]

    def filter_by_tm(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Retain probes with melting temperature (Tm) within the specified range.
        """
        return [r for r in probes if self.tm_min <= mt.Tm_NN(str(r.seq)) <= self.tm_max]

    def filter_by_repeats(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Remove probes with:
        - homopolymers (≥6 nt)
        - dinucleotide repeats (e.g., ATATAT ≥ 3x)
        - trinucleotide repeats (e.g., GCGGCGGCG ≥ 3x)
        """
        def has_repeats(seq: str) -> bool:
            return (
                re.search(r"(A{6,}|T{6,}|G{6,}|C{6,})", seq) or
                re.search(r"((..))\1{2,}", seq) or
                re.search(r"((...))\1{2,}", seq)
            )
        return [r for r in probes if not has_repeats(str(r.seq))]

    def filter_by_structure(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Remove probes with secondary structures more stable than ΔG threshold
        (calculated using RNAfold from ViennaRNA).
        """
        filtered = []
        for record in probes:
            with tempfile.NamedTemporaryFile(mode="w+", suffix=".fa") as temp_in:
                SeqIO.write([record], temp_in, "fasta")
                temp_in.flush()
                try:
                    result = subprocess.run(
                        ["RNAfold", "--noPS", temp_in.name],
                        capture_output=True, text=True, check=True
                    )
                    lines = result.stdout.strip().split("\n")
                    if len(lines) >= 2:
                        dg_match = re.search(r"\(([-\d\.]+)\)", lines[1])
                        if dg_match:
                            dg = float(dg_match.group(1))
                            if dg >= self.dg_threshold:
                                filtered.append(record)
                except Exception as e:
                    print(f"[WARN] RNAfold failed on {record.id}: {e}")
        return filtered

    def apply_all(self) -> None:
        """
        Apply all enabled filters in sequence to input FASTA:
        GC → Tm → Repeats → Structure (optional)
        """
        print(f"[INFO] Reading probes from {self.input_fasta}")
        probes = list(SeqIO.parse(self.input_fasta, "fasta"))
        print(f"[INFO] Total probes loaded: {len(probes)}")

        probes = self.filter_by_gc(probes)
        print(f"[OK] After GC filter: {len(probes)} probes remain")

        probes = self.filter_by_tm(probes)
        print(f"[OK] After Tm filter: {len(probes)} probes remain")

        if not self.allow_repeats:
            probes = self.filter_by_repeats(probes)
            print(f"[✓] After repeat filter: {len(probes)} probes remain")

        if self.structure_filter:
            probes = self.filter_by_structure(probes)
            print(f"[✓] After structure filter: {len(probes)} probes remain")

        SeqIO.write(probes, self.output_fasta, "fasta")
        print(f"[OK] Filtered probes saved to {self.output_fasta}")
