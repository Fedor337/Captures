import os
import urllib.request
import gzip
import shutil
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
            print(f"[\u2713] Уже существует: {destination_path}")
            return

        scheme = url.split("://")[0]
        try:
            if scheme in ("http", "https"):
                print(f"[\u2193] Скачиваем: {url}")
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
                print(f"[\u2193] Скачиваем FTP: {url}")
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

            print(f"[\u2713] Скачано: {destination_path}")
        except Exception as e:
            print(f"[!] Ошибка скачивания: {e}")
            raise

    @staticmethod
    def gunzip_file(input_path: Path, output_path: Path) -> None:
        if output_path.exists():
            print(f"[\u2713] Уже распакован: {output_path}")
            return
        try:
            print(f"[\u21AA] Распаковка: {input_path.name}")
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
            print(f"[\u2713] Распаковано: {output_path}")
        except Exception as e:
            print(f"[!] Ошибка распаковки: {e}")
            raise

    def download_and_extract(self, force_download=False) -> None:
        if not self.gtf.exists() or force_download:
            if not self.gtf_gz.exists() or force_download:
                self.download_file(self.gtf_url, self.gtf_gz)
            self.gunzip_file(self.gtf_gz, self.gtf)
        else:
            print(f"[\u2713] GTF уже распакован: {self.gtf}")
        self.update_progress("GTF загружен и распакован")

        if not self.genome.exists() or force_download:
            if not self.genome_gz.exists() or force_download:
                self.download_file(self.genome_url, self.genome_gz)
            self.gunzip_file(self.genome_gz, self.genome)
        else:
            print(f"[\u2713] Геном уже распакован: {self.genome}")
        self.update_progress("Геном загружен и распакован")

    def get_fasta_chrom_format(self) -> str:
        try:
            with open(self.genome, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        return 'chr' if line[1:].startswith('chr') else ''
        except Exception as e:
            print(f"[!] Ошибка при определении формата хромосом: {e}")
        return ''

    def index_with_bwa(self, force=False) -> None:
        index_files = [self.genome.with_suffix(suffix) for suffix in ['.amb', '.ann', '.bwt', '.pac', '.sa']]
        if all(f.exists() for f in index_files) and not force:
            print(f"[\u2713] Индекс BWA уже существует.")
            self.update_progress("BWA индекс уже существует")
            return
        print(f"[\U0001F527] Строим индекс BWA...")
        try:
            subprocess.run(["bwa", "index", str(self.genome)], check=True)
            print(f"[\u2713] BWA индекс готов.")
            self.update_progress("BWA индекс построен")
        except Exception as e:
            print(f"[!] Ошибка индексации BWA: {e}")
            raise

    def extract_brca_exons(self) -> None:
        print(f"[\U0001F4CD] Извлекаем координаты экзонов BRCA1/2...")
        chrom_prefix = self.get_fasta_chrom_format()

        df = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
        df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "info"]
        exons = df[(df["feature"] == "exon") & (df["info"].str.contains('gene_name \"BRCA1\"|gene_name \"BRCA2\"'))].copy()
        exons["gene"] = exons["info"].str.extract(r'gene_name \"([^\"]+)\"')

        if chrom_prefix:
            exons["chr"] = exons["chr"].astype(str).apply(lambda c: c if c.startswith('chr') else f'chr{c}')
        else:
            exons["chr"] = exons["chr"].astype(str).apply(lambda c: c.replace('chr', ''))

        bed_df = exons[["chr", "start", "end", "gene"]].copy()
        bed_df["start"] = bed_df["start"].astype(int) - 1
        bed_df = bed_df.sort_values(by=["chr", "start"])
        bed_df.to_csv(self.bed, sep='\t', header=False, index=False)
        print(f"[\u2713] Сохранено в BED: {self.bed}")
        self.update_progress("Экзоны BRCA извлечены")

    def extract_sequences_bedtools(self) -> None:
        print(f"[\U0001F9EC] Извлекаем последовательности экзонов через bedtools...")
        cmd = [
            "bedtools", "getfasta",
            "-fi", str(self.genome),
            "-bed", str(self.bed),
            "-fo", str(self.exons_fa),
            "-name"
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"[\u2713] Секвенции сохранены в: {self.exons_fa}")
            self.update_progress("Последовательности экзонов извлечены")
        except Exception as e:
            print(f"[!] Ошибка bedtools getfasta: {e}")
            raise

    def prepare_all(self, force_download=False, force_preparing=False) -> None:
        self.download_and_extract(force_download=force_download)
        self.index_with_bwa(force=force_preparing)
        self.extract_brca_exons()
        self.extract_sequences_bedtools()
        print("[\u2705] Все этапы подготовки завершены")

