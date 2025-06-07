from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.lc_pattern import low_complexity_regions
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
                 dg_threshold: float = -9.0,
                 homopolymer_threshold: int = 5,
                 tandem_min_repeats: int = 3,
                 enable_palindromes: bool = True,
                 enable_low_complexity: bool = True,
                 palindrome_min_length: int = 6):
        self.input_fasta = Path(input_fasta)
        self.output_fasta = Path(output_fasta)
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.allow_repeats = allow_repeats
        self.structure_filter = structure_filter
        self.dg_threshold = dg_threshold

        self.homopolymer_threshold = homopolymer_threshold
        self.tandem_min_repeats = tandem_min_repeats
        self.enable_palindromes = enable_palindromes
        self.enable_low_complexity = enable_low_complexity
        self.palindrome_min_length = palindrome_min_length

    def filter_by_gc(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Retain probes with GC content within the specified range.
        """
        filtered = []

        for record in probes:
            seq = record.seq.upper()
            gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
            if self.gc_min <= gc <= self.gc_max:
                filtered.append(record)

        return filtered

    def filter_by_tm(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Retain probes with melting temperature (Tm) within the specified range.
        """
        return [r for r in probes if self.tm_min <= mt.Tm_NN(r.seq) <= self.tm_max]

    def filter_by_repeats(self, probes: List[SeqRecord]) -> List[SeqRecord]:
        """
        Remove probes with:
        - homopolymers (e.g., AAAAAA)
        - tandem repeats (e.g., ATATAT or GCGGCGGCG)
        - palindromes (e.g., GAATTC)
        - low-complexity regions (via Biopython)
        """
        return [r for r in probes if not self._has_repeats(r.seq)]

    def _has_repeats(self, seq: Seq) -> bool:
        s = str(seq).upper()

        # Homopolymers (e.g., AAAAA)
        homopolymer_re = (r"(A{{{0},}}|T{{{0},}}|G{{{0},}}|C{{{0},}})").format(self.homopolymer_threshold)
        if re.search(homopolymer_re, s):
            return True

        # Tandem repeats (e.g., ATATAT or GCGGCGGCG)
        for size in range(2, 6):
            pattern = re.compile(rf"((\w{{{size}}}))\2{{{self.tandem_min_repeats - 1},}}")
            if pattern.search(s):
                return True

        # Palindromes
        if self.enable_palindromes:
            for i in range(len(seq) - self.palindrome_min_length + 1):  # minimal palindrome length = 6
                for l in range(self.palindrome_min_length, len(seq) - i + 1):
                    fragment = seq[i:i + l]
                    if fragment == fragment.reverse_complement():
                        return True

        # Low complexity (via Biopython DUST-like filter)
        if self.enable_low_complexity:
            if low_complexity_regions(seq):
                return True

        return False

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
            print(f"[OK] After repeat filter: {len(probes)} probes remain")

        if self.structure_filter:
            probes = self.filter_by_structure(probes)
            print(f"[OK] After structure filter: {len(probes)} probes remain")

        SeqIO.write(probes, self.output_fasta, "fasta")
        print(f"[OK] Filtered probes saved to {self.output_fasta}")
