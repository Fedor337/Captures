import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from probe_generator import ProbeGenerator

def test_make_probes_short_sequence():
    """
    Test that a single short exon (< probe length) produces exactly one probe
    with the same sequence.
    """
    pg = ProbeGenerator()
    probes = pg.make_probes("A" * 100, "short_exon")
    assert len(probes) == 1
    assert str(probes[0].seq) == "A" * 100

def test_make_probes_long_sequence():
    """
    Test that a longer exon generates the correct number of overlapping probes
    with the correct length and overlap.
    """
    pg = ProbeGenerator(probe_length=120, max_step=60)
    probes = pg.make_probes("A" * 300, "long_exon")
    # Expected: 4 probes covering 300 nt with step ≤ 60
    assert len(probes) == 4
    for probe in probes:
        assert len(probe.seq) == 120

def test_generate_all(tmp_path):
    """
    Integration test: generates probes from multiple exon sequences,
    writes to a FASTA file, and verifies content.
    """
    # Create mock input FASTA
    fasta_path = tmp_path / "test_exons.fa"
    records = [
        SeqRecord(Seq("A" * 130), id="exon1"),
        SeqRecord(Seq("C" * 250), id="exon2")
    ]
    SeqIO.write(records, fasta_path, "fasta")

    output_path = tmp_path / "test_probes.fa"
    pg = ProbeGenerator(input_fasta=fasta_path, output_fasta=output_path)
    pg.generate_all()

    assert output_path.exists()
    probes = list(SeqIO.parse(output_path, "fasta"))
    assert all(len(p.seq) == 120 or len(p.seq) == 130 for p in probes)
    assert any("exon1" in p.id for p in probes)
    assert any("exon2" in p.id for p in probes)
    assert len(probes) >= 3
