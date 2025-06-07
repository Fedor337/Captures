# BRCA1/2 Probe Design Pipeline

This repository contains a Python-based pipeline for designing oligonucleotide probes targeting the BRCA1 and BRCA2 genes in the human genome (reference: hs37d5). The pipeline automates:
- downloading and preprocessing genome and annotation files,
- extracting exon coordinates,
- generating overlapping probes,
- and filtering probes by GC content, melting temperature (Tm), repeats, and predicted secondary structure (ΔG).

The resulting probe set is suitable for hybridization-based target enrichment in NGS experiments.

---

## 🧬 Features

- Automated download and extraction of genome (FASTA) and annotation (GTF)
- Extraction of BRCA1/2 exon coordinates into BED format
- Sequence retrieval using `bedtools getfasta`
- Generation of overlapping probes (default: 120 nt, step ≤ 60 nt)
- Filtering by:
  - GC content (e.g., 40–60%)
  - Melting temperature (Tm)
  - Repeats (homopolymers, di-/tri-nucleotide patterns, palindromes)
  - Predicted secondary structure (ΔG via RNAfold)
- Executable via Python script (`main.py`) with command-line arguments
- Tested and modular structure for future extension (alignment, reporting)

---

## 📦 Requirements

### 🔹 System Requirements
- Python ≥ 3.8
- [`bedtools`](https://bedtools.readthedocs.io/) ≥ 2.30 — for FASTA extraction
- [`RNAfold`](https://www.tbi.univie.ac.at/RNA/) — optional, used for secondary structure filtering (--structure-filter)

> Install bedtools:
> ```bash
> sudo apt install bedtools        # Debian/Ubuntu
> brew install bedtools            # macOS
> conda install -c bioconda bedtools
> ```

> Install RNAfold (ViennaRNA):
> ```bash
> sudo apt install vienna-rna      # Debian/Ubuntu
> brew install viennarna           # macOS
> conda install -c bioconda viennarna
> ````

### 🔹 Python Packages
Install using pip:

```bash
pip install -r requirements.txt
```

---

## ⚡ Quick Start (command-line)

Run the full pipeline using the provided `main.py` script:

```bash
python main.py
```

#### Optional arguments:

- `--force-download` – re-download reference data
- `--force-prep` – regenerate exon BED and FASTA files
- `--input-fasta path` – override exon input FASTA file (default: `data/brca_exons.fa`)
- `--output-fasta path` – override probe output FASTA file (default: `data/brca_probes.fa`)
- `--probe-length N` – set probe length (default: 120)
- `--max-step N` – set max step between probes (default: 60)

#### Filtering options:
- `--gc-min` <float> – minimum GC content (%), default: 40.0
- `--gc-max` <float> – maximum GC content (%), default: 60.0
- `--tm-min` <float> – minimum Tm (°C), default: 65.0
- `--tm-max` <float> – maximum Tm (°C), default: 72.0
- `--no-repeats` – exclude probes with homopolymers or repeat patterns
- `--structure-filter` – enable RNAfold-based secondary structure filtering
- `--dg-threshold` <float> – minimum acceptable ΔG (kcal/mol), default: -9.0

#### Repeat filter parameters:
- `--homopolymer-threshold <int>` – min. length of homopolymers (default: 6)
- `--tandem-min-repeats <int>` – min. number of motif repeats (default: 3)
- `--disable-palindromes` – disable filtering of palindromic sequences
- `--disable-low-complexity` – disable low-complexity region filtering
- `--palindrome-min-length <int>` – min. length for palindromes (default: 6)

### 🧪 Example Invocations

```bash
# 🚀 Run the full pipeline with defaults (download, extract, generate probes)
python main.py

# 🔄 Redownload genome/annotation and regenerate exon/probe data
python main.py --force-download --force-prep

# 🧬 Use custom probe size and step
python main.py --probe-length 100 --max-step 40

# ⚠️ Test edge case: very sparse probes (no overlap)
python main.py --probe-length 200 --max-step 200

# 📁 Output to a different folder
python main.py --output-fasta results/probes_v1.fa

# 💾 Save to timestamped file (e.g., CI/CD or versioning)
python main.py --output-fasta results/probes_$(date +%Y%m%d).fa

# 🧬 Use a custom exon FASTA file (skip exon extraction)
python main.py --input-fasta data/my_exons.fa \
               --output-fasta results/my_probes.fa

# 🛠 Run probe generation only (assumes reference data exists)
python main.py --force-prep

# 🧪 Use in integration tests or CI pipelines
python main.py --input-fasta data/exons_test.fa \
               --output-fasta data/probes_test.fa

# 🐍 Chain with external filters (example for future use)
python main.py --output-fasta tmp/probes_unfiltered.fa && \
python filter_gc.py --input tmp/probes_unfiltered.fa --output probes_gc_filtered.fa

# 🌡 Filter by melting temperature only
python main.py --tm-min 64 --tm-max 70

# 🚫 Remove probes with repeats
python main.py --no-repeats

# 🧱 Set repeat thresholds explicitly
python main.py --homopolymer-threshold 5 --tandem-min-repeats 4

# 🚫 Disable specific repeat filters
python main.py --disable-palindromes --disable-low-complexity

# 💧 Filter by GC and secondary structure (ΔG ≥ -8.0 kcal/mol)
python main.py --gc-min 42 --gc-max 58 --structure-filter --dg-threshold -8.0

# 👩‍🔬 Quick check with shorter probes (e.g., for microarray simulation)
python main.py --probe-length 80 --max-step 40

# 🔬 Full filtering: GC, Tm, repeats, structure
python main.py \\
  --force-download --force-prep \\
  --probe-length 120 --max-step 60 \\
  --gc-min 40 --gc-max 60 \\
  --tm-min 65 --tm-max 72 \\
  --no-repeats --structure-filter --dg-threshold -9.0
  --homopolymer-threshold 6 --tandem-min-repeats 3 \
  --disable-palindromes --disable-low-complexity \
  --palindrome-min-length 6

```

---

## 🧪 Running Tests

Run all tests using `pytest`:

```bash
pytest
```

> Some tests are skipped automatically if `bedtools` or `RNAfold` is not installed.

---

## 🧬 Using from Python (Alternative)

```python
from reference_preparer import ReferencePreparer
from probe_generator import ProbeGenerator
from probe_filter_pipeline import ProbeFilterPipeline
```

### 📥 Step 1: Prepare reference data
```python
rp = ReferencePreparer(
    genome_url="ftp://...",                          # Optional: override default URLs
    annotation_url="ftp://...",
    output_dir="data"                                # Default: "data"
)
rp.prepare_all(force_download=True, force_preparing=True)  # Download, extract, process. Both arguments are False by default
```

### 🧬 Step 2: Generate overlapping probes from BRCA1/2 exons
```python
pg = ProbeGenerator(
    input_fasta="data/brca_exons.fa",
    output_fasta="data/brca_probes.raw.fa",
    probe_length=120,
    max_step=60
)
pg.generate_all()  # Create overlapping probes from exon sequences
```

### 🧹 Step 3: Apply filters to probes
Each filtering step is available as a separate method and returns a filtered list of SeqRecord objects:
```python
pf = ProbeFilterPipeline(
    input_fasta="data/brca_probes.raw.fa",
    output_fasta="data/brca_probes.filtered.fa",
    gc_min=40,
    gc_max=60,
    tm_min=65,
    tm_max=72,
    allow_repeats=False,
    structure_filter=True,
    dg_threshold=-9.0
    homopolymer_threshold=6,
    tandem_min_repeats=3,
    enable_palindromes=False,
    enable_low_complexity=False,
    palindrome_min_length=6
)
pf.apply_all()
```

### 🛠 Optional: Use filters independently
```python
from Bio import SeqIO

probes = list(SeqIO.parse("data/brca_probes.raw.fa", "fasta"))

filtered_gc = pf.filter_by_gc(probes)
filtered_tm = pf.filter_by_tm(filtered_gc)
filtered_final = pf.filter_by_repeats(filtered_tm)
# Optional: structure filtering (requires RNAfold)
filtered_final = pf.filter_by_structure(filtered_rep)

SeqIO.write(filtered_final, "data/brca_probes.manual.fa", "fasta")
```
This is useful if you want to inspect intermediate results or apply filters interactively in notebooks.

### 🧾 Available Filters

| Method                      | Description                                                                 |
|-----------------------------|-----------------------------------------------------------------------------|
| `filter_by_gc(probes)`      | Keep probes with GC content within `gc_min`–`gc_max` (%)                   |
| `filter_by_tm(probes)`      | Keep probes with melting temperature within `tm_min`–`tm_max` (°C)         |
| `filter_by_repeats(probes)` | Remove probes with homopolymers, tandem repeats, palindromes, low-complexity |
| `filter_by_structure(probes)` | Remove probes with strong secondary structure (ΔG < `dg_threshold`, via RNAfold) |

---

## 🗂 File Structure

```
.
├── data/                         # All intermediate and output files
│   ├── hs37d5.fa                 # Reference genome (unzipped)
│   ├── gencode.v19.annotation.gtf  # Gene annotation (unzipped)
│   ├── brca_exons.bed            # BRCA1/2 exon coordinates
│   ├── brca_exons.fa             # Extracted exon sequences
│   ├── brca_probes.raw.fa        # Raw unfiltered probes
│   ├── brca_probes.fa            # Final filtered probes
│   └── ...                       # Future: alignments, reports
├── reference_preparer.py         # Class for downloading and preprocessing reference data
├── probe_generator.py            # Class for generating overlapping probes
├── probe_filter_pipeline.py      # Class for filtering probes (GC, Tm, repeats, ΔG)
├── main.py                       # Command-line entry point
├── requirements.txt
├── .gitignore
├── test_reference_preparer.py
├── test_probe_generator.py
├── test_probe_filter_pipeline.py
├── test_structure_filter.py
```

---

## 📌 To Do

- Add summary report (number of probes filtered at each step)
- Add support for multi-threaded structure filtering
- Add BLAST/BWA alignment step for specificity checking
- Visualize probe tiling across exons
- Add CLI output in JSON or TSV (optional metadata per probe)
- Add Jupyter Notebook wrapper for exploratory use

---

## 📖 License

MIT License. See `LICENSE` file.

---

## ⚖️ Bash Shell Wrapper (Optional)

For convenience, you may use a simple shell script:

```bash
#!/bin/bash

# Run full pipeline with default parameters
echo "[INFO] Starting BRCA1/2 pipeline"
python main.py \\
    --force-download --force-prep \\
    --probe-length 120 --max-step 60 \\
    --gc-min 40 --gc-max 60 \\
    --tm-min 65 --tm-max 72 \\
    --no-repeats --structure-filter --dg-threshold -9.0
  --homopolymer-threshold 6 --tandem-min-repeats 3 \
  --disable-palindromes --disable-low-complexity \
  --palindrome-min-length 6
```

Save this to `run_pipeline.sh`, then run:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

---


## 🧾 Альтернативная ручная инструкция

```bash
# BRCA1/2 Probe Design Pipeline

# Скачиваем аннотацию генома
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz

# Находим строки относящиеся к BRCA1 и BRCA2
grep -E 'BRCA1|BRCA2' gencode.v19.annotation.gtf | grep 'exon' > brca_exons.gtf

# Достаем координаты
awk 'BEGIN{OFS="\t"} {match($0, /gene_name "([^"]+)"/, a); print $1, $4 - 1, $5, a[1]}' brca_exons.gtf > brca_exons.bed

# Убираем дубли и сортируем
sort -k1,1 -k2,2n brca_exons.bed | uniq > brca_exons.sorted.bed

# Скачиваем последовательность генома
wget https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hs37d5.fa

# Индексируем геном
bwa index hs37d5.fa

# Извлекаем последовательности экзонов
bedtools getfasta -fi hs37d5.fa -bed brca_exons.bed -fo brca_exons.fa -name

# Генерация зондов по экзонам (120 nt с шагом ≤ 60)
# Фильтрация по GC 40–60%
# Фильтрация по температуре 65–72°C
# Удаление зондов с повторами (tandem, low-complexity)
# Удаление зондов со вторичной структурой

# Выравнивание на геном
bwa mem hs37d5.fa probes_final.fa > probes_aligned.sam

# Преобразование в BAM, сортировка, индексация
samtools view -Sb probes_aligned.sam > probes_aligned.bam
samtools sort probes_aligned.bam -o probes_aligned.sorted.bam
samtools index probes_aligned.sorted.bam

# Проверка специфичности (1 зонд = 1 экзон)
samtools view probes_aligned.sorted.bam | cut -f1 | sort | uniq -c | awk '$1 == 1'
```
