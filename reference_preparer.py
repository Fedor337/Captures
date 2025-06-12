import os
import urllib.request
import gzip
import shutil
import subprocess
import pandas as pd
import requests
from pathlib import Path

class ReferencePreparer:
    def __init__(self,
                 genome_url: str = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
                 gtf_url: str = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
                 data_dir: str | Path = "data"):
        self.genome_url = genome_url
        self.gtf_url = gtf_url
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

        self.genome_gz = self.data_dir / "hs37d5.fa.gz"
        self.genome = self.data_dir / "hs37d5.fa"
        self.gtf_gz = self.data_dir / "gencode.v19.annotation.gtf.gz"
        self.gtf = self.data_dir / "gencode.v19.annotation.gtf"
        self.bed = self.data_dir / "brca_exons_sorted.bed"
        self.exons_fa = self.data_dir / "brca_exons.fa"
        self.total_steps = 5
        self.current_step = 0

    def update_progress(self, message: str):
        self.current_step += 1
        percent = int((self.current_step / self.total_steps) * 100)
        bar_length = 40
        filled_length = int(bar_length * percent // 100)
        bar = '=' * filled_length + '-' * (bar_length - filled_length)
        print(f"[{bar}] {percent}% - {message}")

    @staticmethod
    def print_download_bar(downloaded, total, name):
        if total > 0:
            percent = downloaded / total * 100
            done = int(40 * downloaded / total)
            bar = '=' * done + '-' * (40 - done)
            print(f"\r[{bar}] {percent:.0f}% - {name}   ", end='', flush=True)

    @staticmethod
    def download_file(url: str, destination_path: Path, chunk_size: int = 8192, force_download: bool = False) -> None:
        if destination_path.exists() and not force_download:
            print(f"[✓] Уже существует: {destination_path}")
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
                    with open(destination_path, 'wb') as out_file:
                        total = response.headers.get("Content-length")
                        total = int(total) if total is not None else 0
                        downloaded = 0
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

            print(f"[✓] Скачано: {destination_path}")
        except Exception as e:
            print(f"[!] Ошибка скачивания: {e}")
            raise

    @staticmethod
    def gunzip_file(input_path: Path, output_path: Path, force_download: bool = False) -> None:
        if output_path.exists() and not force_download:
            print(f"[✓] Уже распакован: {output_path}")
            return
        try:
            print(f"[↪] Распаковка: {input_path.name}")
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
            print(f"[✓] Распаковано: {output_path}")
        except Exception as e:
            print(f"[!] Ошибка распаковки: {e}")
            raise

    def download_and_extract(self, force_download: bool = False) -> None:
        self.download_file(self.gtf_url, self.gtf_gz, force_download=force_download)
        self.gunzip_file(self.gtf_gz, self.gtf, force_download=force_download)
        self.update_progress("GTF загружен и распакован")

        self.download_file(self.genome_url, self.genome_gz, force_download=force_download)
        self.gunzip_file(self.genome_gz, self.genome, force_download=force_download)
        self.update_progress("Геном загружен и распакован")

    def index_with_bowtie2(self, force: bool = False) -> None:
        index_files = [self.genome.with_suffix(f".fa.{s}.bt2") for s in ['1', '2', '3', '4', 'rev.1', 'rev.2']]
        if all(f.exists() for f in index_files) and not force:
            print(f"[✓] Индекс Bowtie2 уже существует.")
            self.update_progress("Bowtie2 индекс уже существует")
            return
        print(f"[🔧] Строим индекс Bowtie2...")
        try:
            subprocess.run(["bowtie2-build", str(self.genome), str(self.genome)], check=True)
            print(f"[✓] Bowtie2 индекс готов.")
            self.update_progress("Bowtie2 индекс построен")
        except Exception as e:
            print(f"[!] Ошибка индексации Bowtie2: {e}")
            raise

    def extract_brca_exons(self, force_preparing: bool = False) -> None:
        if self.bed.exists() and not force_preparing:
            print(f"[✓] BED уже существует: {self.bed}")
            self.update_progress("BED экзонов уже существует")
            return
        print(f"[📍] Извлекаем координаты экзонов BRCA1/2...")
        df = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
        df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "info"]
        exons = df[(df["feature"] == "exon") & (df["info"].str.contains('gene_name \"BRCA1\"|gene_name \"BRCA2\"'))].copy()
        exons["gene"] = exons["info"].str.extract(r'gene_name \"([^\"]+)\"')
        exons["chr"] = exons["chr"].str.replace("^", "")
        bed_df = exons[["chr", "start", "end", "gene"]].copy()
        bed_df["start"] -= 1
        bed_df = bed_df.sort_values(by=["chr", "start"])
        bed_df.to_csv(self.bed, sep='\t', header=False, index=False)
        print(f"[✓] Сохранено в BED: {self.bed}")
        self.update_progress("Экзоны BRCA извлечены")

    def extract_sequences_bedtools(self, force_preparing: bool = False) -> None:
        if self.exons_fa.exists() and not force_preparing:
            print(f"[✓] FASTA уже существует: {self.exons_fa}")
            self.update_progress("FASTA экзонов уже существует")
            return
        print(f"[🧬] Извлекаем последовательности экзонов через bedtools...")
        cmd = [
            "bedtools", "getfasta",
            "-fi", str(self.genome),
            "-bed", str(self.bed),
            "-fo", str(self.exons_fa),
            "-name"
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"[✓] Секвенции сохранены в: {self.exons_fa}")
            self.update_progress("Последовательности экзонов извлечены")
        except Exception as e:
            print(f"[!] Ошибка bedtools getfasta: {e}")
            raise

    def prepare_all(self, force_download=False, force_preparing=False) -> None:
        self.download_and_extract(force_download=force_download)
        self.index_with_bowtie2(force=force_preparing)
        self.extract_brca_exons(force_preparing=force_preparing)
        self.extract_sequences_bedtools(force_preparing=force_preparing)
        print("[✅] Все этапы подготовки завершены")
