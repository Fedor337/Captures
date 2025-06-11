# Цель
Разработать зонды размерами 120нк к экзонам генов BRCA1/2 для дальнейшего проведения целевого обогащения.

# Теоретическая часть

Целевое обогащение (англ. target enrichment) — это молекулярно-биологический метод, используемый для выделения и обогащения определённых участков генома из общего генетического материала. Такой подход позволяет сфокусироваться на генах или регионах, представляющих интерес, игнорируя ненужные участки, и тем самым экономит ресурсы при секвенировании. Этот метод особенно полезен при секвенировании нового поколения (NGS), так как позволяет значительно сократить объём данных, повышая при этом глубину покрытия целевых регионов и точность анализа.

### Существует несколько основных стратегий целевого обогащения:
#### Гибридизационное обогащение:
Использует биотинилированные олигонуклеотидные зонды, комплементарные интересующим регионам.
Зонды захватывают фрагменты ДНК в растворе (in-solution) или на поверхности (on-array).
Комплексы извлекаются с помощью стрептавидин-модифицированных магнитных шариков.
Пример: Agilent SureSelect, IDT xGen.
#### Амплификационное обогащение:
Использует полимеразную цепную реакцию (ПЦР) с парами праймеров, охватывающими интересующие участки.
Подходит для небольших панелей, но менее масштабируем.
#### CRISPR/Cas-зависимое обогащение:
Основано на направленной разрезке внецелевых участков, оставляя нужные фрагменты.
Менее распространено, но находит применение при работе с длинными фрагментами ДНК.

### Принципы проведения дизайна:
Дизайн зондов — ключевая часть гибридизационного обогащения. Неправильный дизайн может привести к снижению эффективности, специфичности и равномерности покрытия.
Основные этапы и параметры:
Выбор целевых регионов:
Геномные координаты интересующих участков (экзоны, UTR, регуляторные области).
Используются базы данных (RefSeq, Ensembl).

### Разработка зондов:
Зонды представляют собой короткие олигонуклеотиды (обычно 80–120 нуклеотидов).
Каждый зонд должен быть комплементарен целевому участку.
Зонды могут перекрываться для плотного покрытия.
Фильтрация по физико-химическим свойствам:
Длина: обычно 120 нт — обеспечивает достаточную специфичность.
Температура плавления (Tm): зонд должен иметь Tm в узком диапазоне (например, 65–72 °C) для стабильной гибридизации.
GC-содержание: слишком низкое или высокое приводит к нестабильности или неспецифичности.
ΔG вторичной структуры: оценивается энергетическая стабильность нежелательных шпилек и димеров (например, ΔG > –9.0 ккал/моль считается допустимым).
Уникальность: зонд не должен гибридизоваться к внецелевым участкам — проверяется с помощью выравнивания (например, BLAST).
Симуляция гибридизации:
Используются биоинформатические инструменты (например, RNAfold, Primer3, NUPACK), чтобы предсказать вторичную структуру и гибридизацию.
Результаты фильтруются по выбранным порогам.
Производство и тестирование:
Отобранные зонды синтезируются и тестируются in vitro или in silico (для проверки охвата и эффективности).


# Конвеер для дизайна зондов BRCA1/2

Этот репозиторий содержит пайплайн для создания олигонуклеотидных зондов, нацеленных на гены BRCA1 и BRCA2 в геноме человека (референс: hs37d5). Пайлайн позволяет:
- загружать и предварительно обрабатывать файлы геномов и аннотаций,
- извлекать координаты экзонов,
- генерировать перекрывающиеся зонды,
- фильтровать зонды по содержанию GC, температуре плавления (Tm), повторам и предсказанной вторичной структуре (ΔG).

Полученный набор зондов подходит для гибридизационного обогащения мишеней в экспериментах NGS.

---

## 🧬 Возможности

- Автоматическая загрузка и извлечение генома (FASTA) и аннотаций (GTF).
- Извлечение координат экзонов BRCA1/2 в формате BED.
- Получение последовательностей с помощью `bedtools getfasta`.
- Генерация перекрывающихся проб (по умолчанию: 120 нт, шаг ≤ 60 нт).
- Фильтрация по:
  - Содержанию GC (например, 40–60%).
  - Температуре плавления (Tm).
  - Повторам (гомополимеры, ди-/тринуклеотидные паттерны, палиндромы).
  - Предсказанной вторичной структуре (ΔG через RNAfold).
- Запуск через Python-скрипт (`main.py`) с аргументами командной строки.
- Модульная структура для будущего расширения (выравнивание, отчеты).

---

## 📦 Requirements


## 🧬 Возможности

- Автоматическая загрузка и извлечение генома (FASTA) и аннотаций (GTF).
- Извлечение координат экзонов BRCA1/2 в формате BED.
- Получение последовательностей с помощью `bedtools getfasta`.
- Генерация перекрывающихся проб (по умолчанию: 120 нт, шаг ≤ 60 нт).
- Фильтрация по:
  - Содержанию GC (например, 40–60%).
  - Температуре плавления (Tm).
  - Повторам (гомополимеры, ди-/тринуклеотидные паттерны, палиндромы).
  - Предсказанной вторичной структуре (ΔG через RNAfold).
- Модульная структура для будущего расширения (выравнивание, отчеты).



#### Необязательные аргументы:

- `--force-download` – re-download reference data
- `--force-prep` – regenerate exon BED and FASTA files
- `--input-fasta path` – override exon input FASTA file (default: `data/brca_exons.fa`)
- `--output-fasta path` – override probe output FASTA file (default: `data/brca_probes.fa`)
- `--probe-length N` – set probe length (default: 120)
- `--max-step N` – set max step between probes (default: 60)

#### Аргументы фильтрации:
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

### 🧪 Пример запуска

# 🔄 Перезагрузка данных и генерация зондов
python main.py --force-download --force-prep

# 🧬 Использование пользовательских длины и шага зондов
python main.py --probe-length 100 --max-step 40

# 📁 Сохранение в другую папку
python main.py --output-fasta results/probes_v1.fa

# 💾 Сохранение в файл с меткой времени
python main.py --output-fasta results/probes_$(date +%Y%m%d).fa

# 🧬 Использование своего FASTA файла экзонов
python main.py --input-fasta data/my_exons.fa \
               --output-fasta results/my_probes.fa

# 🛠 Генерация зондов без загрузки данных
python main.py --force-prep

# 🌡 Фильтрация только по температуре плавления
python main.py --tm-min 64 --tm-max 70

# 🚫 Удаление проб с повторами
python main.py --no-repeats

# 🧱 Настройка параметров повторов
python main.py --homopolymer-threshold 5 --tandem-min-repeats 4

# 🚫 Отключение отдельных фильтров повторов
python main.py --disable-palindromes --disable-low-complexity

# 💧 Фильтрация по GC и вторичной структуре
python main.py --gc-min 42 --gc-max 58 --structure-filter --dg-threshold -8.0

# 🔬 Полная фильтрация: GC, Tm, повторы, структура
python main.py \\
  --force-download --force-prep \\
  --probe-length 120 --max-step 60 \\
  --gc-min 40 --gc-max 60 \\
  --tm-min 65 --tm-max 72 \\
  --no-repeats --structure-filter --dg-threshold -9.0
  --homopolymer-threshold 6 --tandem-min-repeats 3 \
  --disable-palindromes --disable-low-complexity \
  --palindrome-min-length 6

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

## 📖 License

MIT License. See `LICENSE` file.


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
# Удаление зондов с повторами (гомополимерные, тандемные, низкосложные, паллиндромы)
# Удаление зондов со вторичной структурой

# Выравнивание на геном
bwa mem hs37d5.fa probes_final.fa > probes_aligned.sam

# Преобразование в BAM, сортировка, индексация
samtools view -Sb probes_aligned.sam > probes_aligned.bam
samtools sort probes_aligned.bam -o probes_aligned.sorted.bam
samtools index probes_aligned.sorted.bam

# Вывод специфических качественных ридов в отдельный фаста файл и их количества в терминал
samtools view -q 20 -F 4 probes_aligned.sorted.bam | awk '{print ">"$1"\n"$10}' | tee high_quality_probes.fa | grep "^>" | wc -l


```
