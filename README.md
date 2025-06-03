Скачиваем нотацию генома
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip gencode.v19.annotation.gtf.gz
Находим строки относящиеся к BRCA1 и BRCA2 генам
grep -E 'BRCA1|BRCA2' gencode.v19.annotation.gtf | grep 'exon' > brca_exons.gtf
Достаем столбики с координатами экзонов
awk 'BEGIN{OFS="\t"} {match($0, /gene_name "([^"]+)"/, a); print $1, $4 - 1, $5, a[1]}' brca_exons.gtf > brca_exons.bed
Убираем дубли и сортируем по координатам
sort -k1,1 -k2,2n brca_exons.bed | uniq > brca_exons.sorted.bed
Скачиваем файл с последовательностью генома 
wget https://ilmn-dragen-giab-samples.s3.amazonaws.com/FASTA/hs37d5.fa
Индексируем геном 
bwa index hs37d5.fa
По координатам достаем последовательности из hs37d5.fa и получаем файл с последовательностями всех экзонов
bedtools getfasta -fi hs37d5.fa -bed brca_exons.bed -fo brca_exons.fa -name
Запускаем Probes generation, который должен нам по экзонам создать зонды, полностью покрывающие эти экзоны, размерами 120 нуклеотидов с шагом не более 60 нуклеотидов.
Запускаем GC filtered и находим все зонды имеющие состав GC 40-60%
Запускаем Temperature filtered и извлекаем все зонды с температурой плавления 65-72 градуса
Запускаем Repetions filtered и находим все зонды без длинных тандемных и низкосложных повторов.
Запускаем Structure filtered и убираем зонды со сложными вторичными структурами
