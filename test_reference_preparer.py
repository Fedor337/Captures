import gzip
import pytest
import shutil
from unittest.mock import patch
from reference_preparer import ReferencePreparer

def test_download_http_file(tmp_path):
    url = "https://httpbin.org/robots.txt"
    dest = tmp_path / "robots.txt"

    ReferencePreparer.download_file(url, dest)

    assert dest.exists()
    assert dest.stat().st_size > 0

    content = dest.read_text(encoding='utf-8')
    assert "User-agent" in content

def test_gunzip_file(tmp_path):
    original_text = b"Hello, compressed world!\n"

    input_path = tmp_path / "test.txt.gz"
    output_path = tmp_path / "test.txt"

    with gzip.open(input_path, 'wb') as f:
        f.write(original_text)

    ReferencePreparer.gunzip_file(input_path, output_path)

    assert output_path.exists()
    assert output_path.read_bytes() == original_text

def test_extract_brca_exons(tmp_path):
    gtf_text = """\
        chr17\tHAVANA\texon\t43044295\t43044875\t.\t+\t.\tgene_id "ENSG..."; gene_name "BRCA1";
        chr13\tHAVANA\texon\t32315087\t32316461\t.\t-\t.\tgene_id "ENSG..."; gene_name "BRCA2";
        chr1\tHAVANA\texon\t10000\t10100\t.\t+\t.\tgene_id "ENSG..."; gene_name "NOTCH1";
        """
    gtf_path = tmp_path / "test.gtf"
    gtf_path.write_text(gtf_text)

    rp = ReferencePreparer(data_dir=tmp_path)
    rp.gtf = gtf_path
    rp.bed = tmp_path / "output.bed"

    rp.extract_brca_exons()

    output = rp.bed.read_text().strip().splitlines()
    assert len(output) == 2
    assert output[0].startswith("chr13") or output[1].startswith("chr13")
    assert output[0].endswith("BRCA2") or output[1].endswith("BRCA2")

@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools is not installed")
def test_extract_sequences_bedtools(tmp_path):
    """
    Integration test for extract_sequences_bedtools().

    Requires bedtools to be installed. Creates a minimal genome file (FA),
    a minimal BED file with two regions, and checks that the output .fa 
    file is created and contains expected FASTA headers and sequences.
    """
    ref_fa = tmp_path / "ref.fa"
    bed = tmp_path / "regions.bed"
    out_fa = tmp_path / "out.fa"

    # Пример генома (фейковый)
    ref_fa.write_text(""">chr1
AAAACCCCGGGGTTTTAAAACCCCGGGG
""")

    # Пример координат
    bed.write_text("chr1\t0\t10\ttest1\nchr1\t5\t15\ttest2\n")

    rp = ReferencePreparer(data_dir=tmp_path)
    rp.genome = ref_fa
    rp.bed = bed
    rp.exons_fa = out_fa

    rp.extract_sequences_bedtools()

    # Проверка: out.fa существует и содержит FASTA
    lines = out_fa.read_text().splitlines()
    assert lines[0].startswith(">test1")
    assert set("ACGTN") >= set(lines[1])

@pytest.mark.skipif(shutil.which("bedtools") is not None, reason="bedtools is installed — mocking unnecessary")
def test_extract_sequences_mock(tmp_path):
    """
    Unit test for extract_sequences_bedtools() with subprocess mocked.

    Ensures that the bedtools command is called via subprocess.run()
    without requiring actual bedtools installation.
    """
    rp = ReferencePreparer(data_dir=tmp_path)
    rp.genome = tmp_path / "genome.fa"
    rp.bed = tmp_path / "exons.bed"
    rp.exons_fa = tmp_path / "out.fa"

    with patch("subprocess.run") as mock_run:
        mock_run.return_value = None
        rp.extract_sequences_bedtools()
        mock_run.assert_called_once()
