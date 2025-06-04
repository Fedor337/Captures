# BRCA1/2 Probe Design Pipeline

This repository contains a Python-based pipeline for designing oligonucleotide probes targeting the BRCA1 and BRCA2 genes in the human genome (reference: hs37d5). The pipeline includes downloading genome data, extracting exon coordinates, generating overlapping probes, and filtering them by GC content and melting temperature.

---

## 🧬 Features

- Automatic download and decompression of the genome and annotation
- Extraction of BRCA1/2 exon coordinates from GTF
- Conversion to BED format
- FASTA extraction via `bedtools`
- Modular structure for future probe filtering and alignment

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

---

### 🔹 Python Packages

Install using pip:

```bash
pip install -r requirements.txt
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
│   └── ...                   # Future: probe candidates, filtered sets
├── reference_preparer.py     # Main class for downloading and preprocessing
├── test_loader.py            # Pytest-based tests
├── requirements.txt
└── .gitignore
```

---

## ▶️ How to Run

```python
from reference_preparer import ReferencePreparer

rp = ReferencePreparer()
rp.prepare_all()  # Download, extract, process
```

---

## 🧪 Running Tests

```bash
pytest
```

> Some tests are skipped automatically if `bedtools` is not installed.

---

## 📌 To Do

- Add `ProbeGenerator` class for step 3 (window slicing)
- Add GC/Tm filters
- Add genome alignment step (BLAST or BWA)

---

## 📖 License

MIT License. See `LICENSE` file.

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
# (скрипт Probes generation)

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

# Проверка специфичности (один зонд — один экзон)
samtools view probes_aligned.sorted.bam | cut -f1 | sort | uniq -c | awk '$1 == 1'
```
