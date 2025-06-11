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

    @staticmethod
    def download_file(url: str, destination_path: Path, chunk_size: int = 8192, force_download: bool = False) -> None:
        if os.path.exists(destination_path) and not force_download:
            print(f"[OK] Already exists: {destination_path}")
            return

        scheme = url.split("://")[0]

        try:
            if scheme in ("http", "https"):
                print(f"[GET] Downloading HTTP: {url}")
                with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    with destination_path.open('wb') as f:
                        for chunk in r.iter_content(chunk_size=chunk_size):
                            f.write(chunk)
            elif scheme == "ftp":
                print(f"[GET] Downloading FTP: {url}")
                with urllib.request.urlopen(url) as response, destination_path.open('wb') as out_file:
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        out_file.write(chunk)
            else:
                raise ValueError(f"Unsupported URL scheme: {scheme}")
            print(f"[OK] Downloaded: {destination_path}")
        except Exception as e:
            print(f"[!] Download failed: {e}")
            raise

    @staticmethod
    def gunzip_file(input_path: Path, output_path: Path, force_download: bool = False) -> None:
        if os.path.exists(output_path) and not force_download:
            print(f"[OK] Already extracted: {output_path}")
            return

        try:
            print(f"[GET] Extracting GunZip: {input_path}")
            with gzip.open(input_path, 'rb') as f_in, output_path.open('wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f"[OK] Extracted to: {output_path}")
        except Exception as e:
            print(f"[!] Extraction failed: {e}")
            raise

    def download_and_extract(self, force_download: bool = False) -> None:
        self.download_file(self.gtf_url, self.gtf_gz, force_download=force_download)
        self.gunzip_file(self.gtf_gz, self.gtf, force_download=force_download)

        self.download_file(self.genome_url, self.genome_gz, force_download=force_download)
        self.gunzip_file(self.genome_gz, self.genome, force_download=force_download)

    def index_with_bowtie2(self, force: bool = False) -> None:
        """
        Indexes the genome FASTA file using bowtie2-build if not already indexed.
        Produces files: hs37d5.fa.1.bt2 ... hs37d5.fa.rev.2.bt2
        """
        index_files = [self.genome.with_suffix(f".fa.{s}.bt2") for s in ['1', '2', '3', '4', 'rev.1', 'rev.2']]
        if all(f.exists() for f in index_files) and not force:
            print(f"[OK] Bowtie2 index already exists.")
            return

        print(f"[GET] Building Bowtie2 index...")
        try:
            subprocess.run(["bowtie2-build", str(self.genome), str(self.genome)], check=True)
            print(f"[OK] Bowtie2 index built.")
        except Exception as e:
            print(f"[!] Bowtie2 index build failed: {e}")
            raise

    def extract_brca_exons(self, force_preparing: bool = False) -> None:
        if self.bed.exists() and not force_preparing:
            print(f"[OK] BED file already exists: {self.bed}")
            return

        print(f"[GET] Extracting BRCA exon coordinates...")
        df = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
        df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "info"]

        exons = df[
            (df["feature"] == "exon") &
            (df["info"].str.contains('gene_name "BRCA1"|gene_name "BRCA2"'))
        ].copy()

        exons["gene"] = exons["info"].str.extract(r'gene_name "([^"]+)"')
        exons["chr"] = exons["chr"].str.replace("^", "")  # just in case

        bed_df = exons[["chr", "start", "end", "gene"]].copy()
        bed_df["start"] = bed_df["start"] - 1  # BED format is 0-based
        bed_df = bed_df.sort_values(by=["chr", "start"])
        bed_df.to_csv(self.bed, sep='\t', header=False, index=False)
        print(f"[OK] BED saved: {self.bed}")

    def extract_sequences_bedtools(self, force_preparing: bool = False) -> None:
        if self.exons_fa.exists() and not force_preparing:
            print(f"[OK] FASTA already exists: {self.exons_fa}")
            return

        print(f"[GET] Extracting exon sequences using bedtools...")
        cmd = [
            "bedtools", "getfasta",
            "-fi", str(self.genome),
            "-bed", str(self.bed),
            "-fo", str(self.exons_fa),
            "-name"
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"[OK] Extracted sequences to: {self.exons_fa}")
        except Exception as e:
            print(f"[!] bedtools getfasta failed: {e}")
            raise

    def prepare_all(self, force_download=False, force_preparing=False) -> None:
        self.download_and_extract(force_download=force_download)
        self.index_with_bowtie2(force=force_preparing)  # 🧬 индексация
        self.extract_brca_exons(force_preparing=force_preparing)
        self.extract_sequences_bedtools(force_preparing=force_preparing)
