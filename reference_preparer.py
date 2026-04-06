import os
import urllib.request
import gzip
import subprocess
import pandas as pd
import requests
from pathlib import Path
import math

class ReferencePreparer:
    def __init__(self,
                 genome_url: str = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
                 gtf_url: str = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz",
                 data_dir: str | Path = "data"):
        self.genome_url = genome_url
        self.gtf_url = gtf_url
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

        self.genome_gz = self.data_dir / "hs37d5.fa.gz"
        self.genome = self.data_dir / "hs37d5.fa"
        self.genome_chr = self.data_dir / "hs37d5_chr.fa"
        self.gtf_gz = self.data_dir / "gencode.v43.annotation.gtf.gz"
        self.gtf = self.data_dir / "gencode.v43.annotation.gtf"
        self.bed = self.data_dir / "brca_exons_sorted.bed"
        self.exons_fa = self.data_dir / "brca_exons.fa"
        self.total_steps = 5
        self.current_step = 0

    def update_progress(self, message: str):
        self.current_step += 1
        percent = min((self.current_step / self.total_steps) * 100, 100)
        self.print_bar(percent, message)

    @staticmethod
    def print_bar(percent: float, label: str):
        bar_length = 40
        percent = min(percent, 100)
        filled_length = int(bar_length * percent // 100)
        bar = '=' * filled_length + '-' * (bar_length - filled_length)
        percent_str = f"{percent:3.0f}%"
        print(f"\r[{bar}] {percent_str} - {label:<35}", flush=True)

    @staticmethod
    def print_download_bar(downloaded, total, name):
        if total > 0:
            percent = min(downloaded / total * 100, 100)
            done = min(int(40 * downloaded / total), 40)
            bar = '=' * done + '-' * (40 - done)
            percent_str = f"{percent:3.0f}%"
            print(f"\r[{bar}] {percent_str} - {name:<35}", end='', flush=True)

    @staticmethod
    def download_file(url: str, destination_path: Path, chunk_size: int = 8192) -> None:
        if destination_path.exists():
            print(f"Уже существует: {destination_path}")
            return

        scheme = url.split("://")[0]
        try:
            if scheme in ("http", "https"):
                print(f"[↓] Скачиваем: {url}")
                with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    total = int(r.headers.get('Content-Length', 0))
                    downloaded = 0
                    with open(destination_path, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=chunk_size):
                            f.write(chunk)
                            downloaded += len(chunk)
                            ReferencePreparer.print_download_bar(downloaded, total, destination_path.name)
                print()
            elif scheme == "ftp":
                print(f"[↓] Скачиваем FTP: {url}")
                with urllib.request.urlopen(url) as response:
                    meta = response.info()
                    total = int(meta.get("Content-Length", 0))
                    downloaded = 0
                    with open(destination_path, 'wb') as out_file:
                        while True:
                            chunk = response.read(chunk_size)
                            if not chunk:
                                break
                            out_file.write(chunk)
                            downloaded += len(chunk)
                            ReferencePreparer.print_download_bar(downloaded, total, destination_path.name)
                print()
            else:
                raise ValueError(f"Unsupported URL scheme: {scheme}")
            print(f"Скачано: {destination_path}")
        except Exception as e:
            print(f"Ошибка скачивания: {e}")
            raise

    @staticmethod
    def gunzip_file(input_path: Path, output_path: Path) -> None:
        if output_path.exists():
            print(f"Уже распакован: {output_path}")
            return
        try:
            print(f"Распаковка: {input_path.name}")
            total = os.path.getsize(input_path)
            processed = 0
            with gzip.open(input_path, 'rb') as f_in, open(output_path, 'wb') as f_out:
                while True:
                    chunk = f_in.read(8192)
                    if not chunk:
                        break
                    f_out.write(chunk)
                    processed += len(chunk)
                    ReferencePreparer.print_download_bar(processed, total, output_path.name)
            print()
            print(f"Распаковано: {output_path}")
        except Exception as e:
            print(f"Ошибка распаковки: {e}")
            raise

    def download_and_extract(self, force_download=False) -> None:
        if not self.gtf.exists() or force_download:
            if not self.gtf_gz.exists() or force_download:
                self.download_file(self.gtf_url, self.gtf_gz)
            self.gunzip_file(self.gtf_gz, self.gtf)
        else:
            print(f"GTF уже распакован: {self.gtf}")
        self.update_progress("GTF загружен и распакован")

        if not self.genome.exists() or force_download:
            if not self.genome_gz.exists() or force_download:
                self.download_file(self.genome_url, self.genome_gz)
            self.gunzip_file(self.genome_gz, self.genome)
        else:
            print(f"Геном уже распакован: {self.genome}")
        self.update_progress("Геном загружен и распакован")

        # Убедимся, что hs37d5_chr.fa с префиксами chr создан
        if not self.genome_chr.exists():
            try:
                print("Добавляем префикс 'chr' к заголовкам FASTA...")
                with open(self.genome, 'r') as f_in, open(self.genome_chr, 'w') as f_out:
                    for line in f_in:
                        if line.startswith('>') and not line.startswith('>chr'):
                            f_out.write('>chr' + line[1:])
                        else:
                            f_out.write(line)
                print(f"[✓] Создан файл: {self.genome_chr}")
            except Exception as e:
                print(f"[!] Ошибка преобразования chr: {e}")
                raise

    def get_fasta_chrom_format(self) -> str:
        try:
            with open(self.genome_chr if self.genome_chr.exists() else self.genome, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        return 'chr' if line[1:].startswith('chr') else ''
        except Exception as e:
            print(f"Ошибка при определении формата хромосом: {e}")
        return ''

    def index_with_bwa(self, force=False) -> None:
        index_base = self.genome_chr if self.genome_chr.exists() else self.genome
        index_files = [index_base.with_suffix(ext) for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        if all(f.exists() for f in index_files) and not force:
            print(f"Индексы BWA уже существуют для {index_base}")
            self.update_progress("BWA индекс уже существует")
            return
        print(f"Строим индекс BWA для {index_base}...")
        try:
            subprocess.run(["bwa", "index", str(index_base)], check=True)
            print(f"BWA индекс готов.")
            self.update_progress("BWA индекс построен")
        except Exception as e:
            print(f"Ошибка индексации BWA: {e}")
            raise

    def extract_brca_exons(self) -> None:
        print(f"Извлекаем координаты экзонов BRCA1/2...")
        chrom_prefix = 'chr'

        df = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
        df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "info"]
        exons = df[(df["feature"] == "exon") & df["info"].str.contains('gene_name "BRCA1"|gene_name "BRCA2"')].copy()
        exons["gene"] = exons["info"].str.extract(r'gene_name "([^"]+)"')

        exons["chr"] = exons["chr"].astype(str).apply(lambda c: f'chr{c}' if not c.startswith('chr') else c)

        bed_df = exons[["chr", "start", "end", "gene"]].copy()
        bed_df["start"] = bed_df["start"].astype(int) - 1
        bed_df = bed_df.sort_values(["chr","start"])
        bed_df.to_csv(self.bed, sep='\t', header=False, index=False)
        print(f"Сохранено в BED: {self.bed}")
        self.update_progress("Экзоны BRCA извлечены")

    def extract_sequences_bedtools(self) -> None:
        print(f"Извлекаем последовательности экзонов через bedtools...")
        try:
            subprocess.run([
                "bedtools", "getfasta",
                "-fi", str(self.genome_chr),
                "-bed", str(self.bed),
                "-fo", str(self.exons_fa),
                "-name"
            ], check=True)
            print(f"Секвенции сохранены в: {self.exons_fa}")
        except Exception as e:
            print(f"Ошибка bedtools getfasta: {e}")
            raise

        # Удаляем дублирующиеся последовательности
        try:
            seen = set()
            unique = []
            with open(self.exons_fa) as f:
                name, seq = None, []
                for line in f:
                    if line.startswith('>'):
                        if name and ''.join(seq) not in seen:
                            seen.add(''.join(seq)); unique.append((name,''.join(seq)))
                        name = line.strip()
                        seq=[]
                    else:
                        seq.append(line.strip())
                if name and ''.join(seq) not in seen:
                    unique.append((name,''.join(seq)))
            with open(self.exons_fa,'w') as f:
                for nm,sq in unique:
                    f.write(f"{nm}\n{sq}\n")
            print(f"Уникальные экзоны: {len(unique)}")
            self.update_progress(f"Уникальные экзоны ({len(unique)}) сохранены")
        except Exception as e:
            print(f"Ошибка уникализации: {e}")
            raise

    def prepare_all(self, force_download=False, force_preparing=False) -> None:
        self.download_and_extract(force_download=force_download)
        self.index_with_bwa(force=force_preparing)
        self.extract_brca_exons()
        self.extract_sequences_bedtools()
        print("[✅] Все этапы подготовки завершены")


