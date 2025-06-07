import shutil
import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from probe_filter_pipeline import ProbeFilterPipeline

RNAFOLD_EXISTS = shutil.which("RNAfold") is not None

def test_gc_filter_passes_sequence(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("ATGC" * 30), id="probe1")  # GC = 50%
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 1
    assert records[0].id == "probe1"

def test_gc_filter_blocks_low_gc(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("ATATATATATATATATATAT"), id="lowgc")  # GC = 0%
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0

def test_gc_filter_mixed(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    records = [
        SeqRecord(Seq("ATGC" * 30), id="good"),
        SeqRecord(Seq("ATATATATAT"), id="low"),
        SeqRecord(Seq("GCGCGCGCGC"), id="high")
    ]
    SeqIO.write(records, input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path, gc_min=30, gc_max=70)
    pipeline.apply_all()
    result_ids = [rec.id for rec in SeqIO.parse(output_path, "fasta")]
    assert result_ids == ["good"]

def test_tm_filter_keeps_within_range(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("ATGC" * 30), id="ok_tm")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path,
                                   gc_min=0, gc_max=100, tm_min=65, tm_max=72)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 1
    assert records[0].id == "ok_tm"

def test_tm_filter_blocks_outside_range(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("ATATATAT"), id="low_tm")  # short = low Tm
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path,
                                   gc_min=0, gc_max=100, tm_min=65, tm_max=72)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0

def test_repeat_filter_blocks_homopolymer(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("A" * 40), id="homo_A")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path,
                                   gc_min=0, gc_max=100, tm_min=0, tm_max=100)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0

def test_repeat_filter_blocks_dinucleotide(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("ATATATATATATAT"), id="dinuc")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path,
                                   gc_min=0, gc_max=100, tm_min=0, tm_max=100)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0

def test_repeat_filter_blocks_trinucleotide(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("GCAGCAGCAGCA"), id="trinuc")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path, output_fasta=output_path,
                                   gc_min=0, gc_max=100, tm_min=0, tm_max=100)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0

def test_repeat_filter_can_be_disabled(tmp_path):
    input_path = tmp_path / "input.fa"
    output_path = tmp_path / "output.fa"
    record = SeqRecord(Seq("A" * 40), id="allowed_homo")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(input_fasta=input_path,
                                   output_fasta=output_path,
                                   gc_min=0, gc_max=100,
                                   tm_min=0, tm_max=100,
                                   allow_repeats=True)
    pipeline.apply_all()
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 1
    assert records[0].id == "allowed_homo"

@pytest.mark.skipif(not RNAFOLD_EXISTS, reason="RNAfold not installed")
def test_structure_filter_excludes_strong_hairpin(tmp_path):
    input_path = tmp_path / "hairpin.fa"
    output_path = tmp_path / "filtered.fa"

    # Sequence likely to form a stable hairpin (low ΔG)
    seq = Seq("GCGCGCGCGCGCGCGCGCGC")
    record = SeqRecord(seq, id="hairpin")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(
        input_fasta=input_path,
        output_fasta=output_path,
        gc_min=0, gc_max=100,
        tm_min=0, tm_max=100,
        allow_repeats=True,
        structure_filter=True,
        dg_threshold=-5.0
    )
    pipeline.apply_all()

    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0

@pytest.mark.skipif(not RNAFOLD_EXISTS, reason="RNAfold not installed")
def test_structure_filter_allows_weak_structure(tmp_path):
    input_path = tmp_path / "flat.fa"
    output_path = tmp_path / "filtered.fa"

    # Sequence with low structure potential (random AT-rich)
    seq = Seq("ATATATATATATATATATAT")
    record = SeqRecord(seq, id="flat")
    SeqIO.write([record], input_path, "fasta")

    pipeline = ProbeFilterPipeline(
        input_fasta=input_path,
        output_fasta=output_path,
        gc_min=0, gc_max=100,
        tm_min=0, tm_max=100,
        allow_repeats=True,
        structure_filter=True,
        dg_threshold=-5.0
    )
    pipeline.apply_all()

    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 1
    assert records[0].id == "flat"
