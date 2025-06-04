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

---

# BRCA1/2 Probe Design Pipeline

- Скачиваем нотацию генома
- wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
- gunzip gencode.v19.annotation.gtf.gz
- Находим строки относящиеся к BRCA1 и BRCA2 генам
- grep -E 'BRCA1|BRCA2' gencode.v19.annotation.gtf | grep 'exon' > brca_exons.gtf
- Достаем столбики с координатами экзонов
- awk 'BEGIN{OFS="\t"} {match($0, /gene_name "([^"]+)"/, a); print $1, $4 - 1, $5, a[1]}' brca_exons.gtf > brca_exons.bed
- Убираем дубли и сортируем по координатам
- sort -k1,1 -k2,2n brca_exons.bed | uniq > brca_exons.sorted.bed
- Скачиваем файл с последовательностью генома 
- wget https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hs37d5.fa
- Индексируем геном 
- bwa index hs37d5.fa
- По координатам достаем последовательности из hs37d5.fa и получаем файл с последовательностями всех экзонов
- bedtools getfasta -fi hs37d5.fa -bed brca_exons.bed -fo brca_exons.fa -name
- Запускаем Probes generation, который должен нам по экзонам создать зонды, полностью покрывающие эти экзоны, размерами 120 нуклеотидов с шагом не более 60 нуклеотидов.
- Запускаем GC filtered и находим все зонды имеющие состав GC 40-60%
- Запускаем Temperature filtered и извлекаем все зонды с температурой плавления 65-72 градуса
- Запускаем Repetions filtered и находим все зонды без длинных тандемных и низкосложных повторов.
- Запускаем Structure filtered и убираем зонды со сложными вторичными структурами
- Выравниваем зонды на геном
- bwa mem hs37d5.fa probes_final.fa > probes_aligned.sam
- Преобразкем sam в bam и сортируем
- samtools view -Sb probes_aligned.sam > probes_aligned.bam
- samtools sort probes_aligned.bam -o probes_aligned.sorted.bam
- samtools index probes_aligned.sorted.bam
- Проверяем специфичность(один зонд-один экзон)
- samtools view probes_aligned.sorted.bam | cut -f1 | sort | uniq -c | awk '$1 == 1'

---
