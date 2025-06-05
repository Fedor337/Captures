# BRCA1/2 Probe Design Pipeline

This repository contains a Python-based pipeline for designing oligonucleotide probes targeting the BRCA1 and BRCA2 genes in the human genome (reference: hs37d5). The pipeline includes downloading genome data, extracting exon coordinates, generating overlapping probes, and filtering them by GC content and melting temperature.

---

## 🧬 Features

- Automated download and extraction of genome (FASTA) and annotation (GTF)
- Extraction of BRCA1/2 exon coordinates into BED format
- Sequence retrieval using `bedtools getfasta`
- Generation of overlapping probes (default: 120 nt, step ≤ 60 nt)
- Executable via Python script (`main.py`) with command-line arguments
- Tested and modular structure for future extension (filtering, alignment)

---

## 📦 Requirements

### 🔹 System Requirements
- Python ≥ 3.8
- [`bedtools`](https://bedtools.readthedocs.io/) ≥ 2.30

> Install bedtools with one of the following:
>
> ```bash
> sudo apt install bedtools        # Debian/Ubuntu
> brew install bedtools            # macOS
> conda install -c bioconda bedtools
> ```

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


### Example Invocations

```bash
# 🚀 Run the full pipeline with defaults (download, extract, generate probes)
python main.py

# 🔄 Redownload genome/annotation and regenerate exon/probe data
python main.py --force-download --force-prep

# 🔬 Change probe length and step size
python main.py --probe-length 100 --max-step 50 \
               --output-fasta data/probes_len100_step50.fa

# 📁 Output to a different folder
python main.py --output-fasta results/probes_v1.fa

# 🧬 Use a custom exon FASTA file (skip exon extraction)
python main.py --input-fasta data/my_exons.fa \
               --output-fasta results/my_probes.fa

# ⚠️ Test edge case: very sparse probes (no overlap)
python main.py --probe-length 200 --max-step 200

# 🧪 Use in integration tests / pipelines
python main.py --input-fasta data/exons_test.fa \
               --output-fasta data/probes_test.fa

# 🛠 Run probe generation only, skip all downloads (assumes data exists)
python main.py --force-prep

# 🐍 Chain with additional filters (future): generate probes, pass to next step
python main.py --output-fasta tmp/probes_unfiltered.fa && \
python filter_gc.py --input tmp/probes_unfiltered.fa --output probes_gc_filtered.fa

# 👩‍🔬 Quick check with shorter probes (e.g., for tiling microarray simulation)
python main.py --probe-length 80 --max-step 40

# 💾 Save to timestamped file (e.g., CI/CD or versioning)
python main.py --output-fasta results/probes_$(date +%Y%m%d).fa
```

---

## 🧪 Running Tests

Run all tests using `pytest`:

```bash
pytest
```

> Some tests are skipped automatically if `bedtools` is not installed.

---

## 🧬 Using from Python (Alternative)

```python
from reference_preparer import ReferencePreparer
from probe_generator import ProbeGenerator
```

### Step 1: Prepare reference data
```python
rp = ReferencePreparer()
rp.prepare_all()  # Download, extract, process
```

### Step 2: Generate probes from BRCA1/2 exons
```python
pg = ProbeGenerator()
pg.generate_all()  # Create overlapping probes from exon sequences
```

---

## 🗂 File Structure

```
.
├── data/                     # All generated data and intermediate files
│   ├── hs37d5.fa             # Reference genome (unzipped)
│   ├── gencode.v19.annotation.gtf  # Gene annotations (unzipped)
│   ├── brca_exons.bed        # BRCA1/2 exon coordinates
│   ├── brca_exons.fa         # Extracted exon sequences
│   ├── brca_probes.fa        # Designed overlapping probes
│   └── ...                   # Future: filtered sets, alignment output
├── reference_preparer.py     # Main class for downloading and preprocessing
├── probe_generator.py        # Class for generating overlapping probes
├── test_reference_preparer.py
├── test_probe_generator.py
├── main.py
├── requirements.txt
└── .gitignore
```

---

## 📌 To Do

- Add GC/Tm/repeats/structure filtering
- Add genome alignment step (BLAST or BWA)
- Add CLI or Jupyter runner

---

## 📖 License

MIT License. See `LICENSE` file.

---

## ⚖️ Bash Shell Wrapper

For convenience, you may use a simple shell script:

```bash
#!/bin/bash

# Run full pipeline with default parameters
echo "[INFO] Starting BRCA1/2 pipeline"
python main.py --force-download --force-prep \
               --probe-length 120 --max-step 60 \
               --output-fasta data/brca_probes.fa
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
