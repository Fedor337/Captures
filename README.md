# Цель
Разработать зонды размерами 120нк к экзонам генов BRCA1/2 для дальнейшего проведения целевого обогащения.

# Теоретическая часть

Целевое обогащение — это молекулярно-биологический метод, используемый для избирательной концентрации определённых фрагментов ДНК или РНК перед секвенированием. Основной задачей метода является выделение участков, представляющих интерес, с помощью гибридизационных зондов или мультиплексной амплификации, исключая остальной генетический материал. Это позволяет существенно снизить объём обрабатываемых данных, увеличить глубину покрытия выбранных регионов, повысить чувствительность обнаружения мутаций, а также существенно снизить затраты на реагентах. Метод широко применяется в диагностике, экзомном и таргетном секвенировании, а также в исследованиях генетической вариабельности и экспрессии.

### Существует несколько основных стратегий целевого обогащения:
#### Гибридизационное обогащение:
Использует биотинилированные олигонуклеотидные зонды, комплементарные интересующим регионам.
Зонды захватывают фрагменты ДНК в растворе (in-solution) или на поверхности (on-array).
Комплексы извлекаются с помощью стрептавидин-модифицированных магнитных шариков.
Пример: Agilent SureSelect, IDT xGen.
Преимущества: Более равномерное покрытие, высокая специфичность, возможность работать с большим числом целей.
#### Амплификационное обогащение:
Использует полимеразную цепную реакцию (ПЦР) с парами праймеров, охватывающими интересующие участки.
Преимущества: Быстрота, низкая стоимость, подходит для деградированной ДНК.
Недостатки: Ограничения по числу мишеней, неравномерное покрытие;
#### CRISPR/Cas-зависимое обогащение:
Основано на направленной разрезке внецелевых участков, оставляя нужные фрагменты.
Менее распространено, но находит применение при работе с длинными фрагментами ДНК.

### Принципы проведения дизайна:
Дизайн зондов — ключевая часть гибридизационного обогащения. Неправильный дизайн может привести к снижению эффективности, специфичности и равномерности покрытия.
1. Определение цели и выбор референсной последовательности.
2. Выбор стратегии захвата, параметров дизайна, программ и инструментов.
3. Фильтрация и валидация зондов in silico.
4. Генерация списка зондов с координатами, последовательностями, характеристиками.
5. Выбор формата заказа и поставщика зондов.
6. Оптимизация и лабораторное тестирование.

### Разработка зондов и ход работы:
1. Выбираем цель обогащения, в нашем случае это экзоны генов BRCA1/2 гена человека. Для это мы используем геномную последовательность человека hs37d5, скачаную из хранилища, принадлежащее компании Illumina и нотацию, которую взяли с сайта EBI.
2. С помощью команд терминала отбираем геномные координаты интересующих нас участков, а именно экзоны, и записываем их в отдельный .gtf файл.
3. Индексиркем наш геном человека с помощью утилиты BWA.
4. По координатам достаем из генома последовательности экзонов с помощью утилиты bedtools.
5. Шаги 1-4 можно пропустить запустив скрипт main.py, который проведет скачивание, индексирование и извлечение экзонных последовательностей автоматически.
6. Разрезаем экзоны на отдельные олигонуклеотидные последовательности размерами 120 нуклеотидов перекрытием от 50% (здесь будет представлено перекрытие на 99%), из которых мы будем отбирать зонды и удаляем дубликаты последовательностей.
7. Дальше проводим фильтрацию:
- Температура плавления (Tm): зонд должен иметь Tm в узком диапазоне, приблизительно 65–72 °C, для стабильной гибридизации (температура является самым важным фактором отбора). Для фильтрации запусаем код temperature_filter.py.
- GC-содержание: ниже 40% или выше 60% приводит к нестабильности или неспецифичности. При создании дизайна, выход при данных параметрах был нулевым, нижнию границу необходимо было опустить до 30% верхнюю до 40%. Запускаем код gc_filter.py.
- Повторы: зонды не должны иметь длинных гомополимерных, тандемных и низкосложных повторов, так как они приводят к низкой специфичности, образованию шпилек и слабой способностью к связыванию. Отбор производился по меньше 5 гомополимерных повторов; энтропии Шеннона 1.8; рамерам мотива для тандема до 5 гомополимеров, повторам для тандема от 3; палиндромам (все параметры идут по умолчанию, но их также можно запускать поотдельности). Запускаем код repetition_filter.py
7. Уникальность: зонд не должен гибридизоваться к внецелевым участкам — проверяется с помощью выравнивания.
8. Производство и тестирование:
Отобранные зонды синтезируются и тестируются in vitro или in silico для проверки охвата и качества.

# Как пользоваться конвейером для дизайна зондов 

## Описание
Этот репозиторий содержит конвейер для создания олигонуклеотидных зондов. Конвейер позволяет:
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
- Запуск через Python-скрипт (`main.py`).
- Модульная структура для будущего расширения (выравнивание, отчеты).

---

## 🗂 Файловая структура

```
.
├── data/                         # Все промежуточные и входные файлы
│   ├── hs37d5.fa                 # Разархивированный референсный геном
│   ├── gencode.v19.annotation.gtf  # Разархивированная аннотация генома
│   ├── brca_exons.bed            # Координаты экзонов BRCA1/2 
│   ├── brca_exons.fa             # Координаты и последовательности экзонов BRCA1/2
│   └── ...                       # В будущем: выравнивания, отчеты
├── reference_preparer.py         # Код для загрузки и предварительной обработки референсных данных
├── probe_generator.py            # Код для генерации перекрывающихся зондов
├── main.py                       # Вход в коммандную строку
├── requirements.txt
├── .gitignore
├── reference_preparer.py
├── temperature_filter.py
├── gc_filter.py
├── dedublicator.py
├── repetitions_filter.py
```

---
## Пример запуска
```
# Запускаем скрипт, который подготавливает дирректории и загружает необходимые для работы последовательности и аннотации
python3 main.py --force-download --force-prep

```
## 📖 License

MIT License. Смотреть `LICENSE` файл.


## 🧾 Альтернативная ручная инструкция

```bash
# Конвейер для генерации и фильтрации зондов к генам

# Скачиваем аннотацию генома
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz 
gunzip gencode.v43.annotation.gtf.gz

# Находим строки относящиеся к BRCA1 и BRCA2
grep -E 'BRCA1|BRCA2' gencode.v43.annotation.gtf | grep 'exon' > brca_exons.gtf

# Достаем координаты
awk 'BEGIN{OFS="\t"} {match($0, /gene_name "([^"]+)"/, a); print $1, $4 - 1, $5, a[1]}' brca_exons.gtf > brca_exons.bed

# Убираем дубли и сортируем
sort -k1,1 -k2,2n brca_exons.bed | uniq > brca_exons_sorted.bed

# Скачиваем последовательность генома
wget https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hs37d5.fa

#Добавляем префиксы chr перед хромосомами в файле hs37d5.fa
sed 's/^>/>chr/' hs37d5.fa > hs37d5_chr.fa

# Индексируем геном
bwa index hs37d5_chr.fa

# Извлекаем последовательности экзонов
bedtools getfasta -fi hs37d5_chr.fa -bed brca_exons.bed -fo brca_exons.fa -name

# Генерация зондов по экзонам (по умолчанию 120 нк с шагом = 1)
python3 probe_generator.py brca_exons.fa brca_probes.fa --probe-length 120 --step 1
# Удаляем дубликаты зондов по последовательностям
python3 dedublicator.py brca_probes.fa brca_probes_unique.fa
# Фильтрация по температуре (по умолчанию 65–72°C)
python3 tm_filter.py brca_probes_unique.fa probes_tm_filtered.fa --tm_min 65 --tm_max 72
# Фильтрация по GC (по умолчанию 40–60%)
python3 gc_filter.py probes_tm_filtered.fa probes_tm_gc_filtered.fa --gc_min 30 --gc_max 40
# Удаление зондов с повторами (гомополимерные, тандемные, низкосложные, паллиндромы)
python3 gc_filter.py  probes_tm_gc_filtered.fa  probes_tm_gc_repetitions_filtered.fa --max-homopolymer 5 --min-entropy 1.8 --tandem-min-repeats 3 --tandem-max-motif 5

# Выравнивание на геном
bwa mem hs37d5_chr.fa probes_tm_gc_repetitions_structure_filtered.fa > probes_aligned.sam

# Преобразование в BAM, сортировка, индексация
samtools view -Sb probes_aligned.sam > probes_aligned.bam
samtools sort probes_aligned.bam -o probes_aligned_sorted.bam
samtools index probes_aligned_sorted.bam

# Вывод специфических качественных ридов в отдельный фаста файл и их количества в терминал
samtools view -q 20 -F 4 probes_aligned_sorted.bam | awk '{print ">"$1"\n"$10}' | tee high_quality_probes.fa | grep "^>" | wc -l

```
