import os
import urllib.request
import gzip
import shutil
import subprocess
import pandas as pd
import requests
from pathlib import Path
# from tqdm import tqdm  # Установи через pip install tqdm


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
        """
        Downloads a file from a given URL (http, https, or ftp) to a local path.

        Supports streaming download for large files and avoids loading entire content into memory.

        Args:
            url (str): The URL to download from. Supports http(s) and ftp.
            dest_path (Path): Path to save the downloaded file.
            chunk_size (int): Number of bytes to read per iteration. Default is 8192 (8 KB).

        Raises:
            ValueError: If the URL scheme is unsupported.
            Exception: If the download fails for any reason (network error, etc.).
        """

        if os.path.exists(destination_path) and not force_download:
            print(f"[OK] Already exists: {destination_path}")
            return

        scheme = url.split("://")[0]

        if scheme in ("http", "https"):
            try:
                with requests.get(url, stream=True) as r:
                    r.raise_for_status()
                    print(f"[GET] Downloading HTTP: {url}")
                    with destination_path.open('wb') as f:
                        for chunk in r.iter_content(chunk_size=chunk_size):
                            f.write(chunk)
                print(f"[OK] Downloaded HTTP: {destination_path}")
            except Exception as e:
                print(f"[!] HTTP download failed: {e}")
                raise

        elif scheme == "ftp":
            try:
                print(f"[GET] Downloading FTP: {url}")
                with urllib.request.urlopen(url) as response, destination_path.open('wb') as out_file:
                                while True:
                                    chunk = response.read(chunk_size)
                                    if not chunk:
                                        break
                                    out_file.write(chunk)
                                print(f"[OK] Downloaded FTP: {destination_path}")
            except Exception as e:
                print(f"[!] FTP download failed: {e}")
                raise

        else:
            raise ValueError(f"Unsupported URL scheme: {scheme}")

    @staticmethod
    def gunzip_file(input_path: Path, output_path: Path, force_download: bool = False) -> None:
        """
        Extracts a .gz file to a specified output file.

        Uses a streaming approach to minimize memory usage. Skips extraction
        if the output file already exists.

        Args:
            input_path (Path): Path to the input .gz file.
            output_path (Path): Path to write the uncompressed output file.

        Raises:
            Exception: If decompression fails (e.g., corrupt .gz file).
        """

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
        """
        Downloads and extracts both the genome FASTA file and GTF annotation file.

        Skips files if already present in the working directory.
        """

        self.download_file(self.gtf_url, self.gtf_gz, force_download=force_download)
        self.gunzip_file(self.gtf_gz, self.gtf, force_download=force_download)

        self.download_file(self.genome_url, self.genome_gz, force_download=force_download)
        self.gunzip_file(self.genome_gz, self.genome, force_download=force_download)

    def extract_brca_exons(self, force_preparing: bool = False) -> None:
        """
        Extracts exon coordinates for BRCA1 and BRCA2 genes from a GTF annotation file,
        and saves them as a BED-format file with columns: chr, start, end, gene.

        This method:
        - reads the GTF file,
        - filters only 'exon' features,
        - selects only rows with gene_name BRCA1 or BRCA2,
        - adjusts start coordinate to 0-based for BED format,
        - and writes the result to a sorted .bed file.

        Output file is saved to: self.bed
        """
        if self.bed.exists() and not force_preparing:
            print(f"[OK] Already exists: {self.bed}")
            return

        print(f"[GET] Extracting BRCA1/2 exon coordinates from: {self.gtf}")
        df = pd.read_csv(self.gtf, sep='\t', comment='#', header=None)
        df.columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "info"]

        exons = df[
            (df["feature"] == "exon") &
            (df["info"].str.contains('gene_name "BRCA1"|gene_name "BRCA2"'))
        ].copy()

        exons["gene"] = exons["info"].str.extract(r'gene_name "([^"]+)"')
        bed_df = exons[["chr", "start", "end", "gene"]].copy()
        bed_df["start"] = bed_df["start"] - 1
        bed_df = bed_df.sort_values(by=["chr", "start"])
        bed_df.to_csv(self.bed, sep='\t', header=False, index=False)
        print(f"[OK] Saved BED: {self.bed}")

    def extract_sequences_bedtools(self, force_preparing: bool = False) -> None:
        if self.exons_fa.exists() and not force_preparing:
            print(f"[OK] Already exists: {self.exons_fa}")
            return

        print(f"[GET] Extracting exon sequences from genome")
        cmd = [
            "bedtools", "getfasta",
            "-fi", self.genome,
            "-bed", self.bed,
            "-fo", self.exons_fa,
            "-name"
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"[OK] Extracted sequences to: {self.exons_fa}")
        except Exception as e:
            print(f"[!] bedtools failed: {e}")
            raise

    def prepare_all(self, force_download=False, force_preparing=False) -> None:
        self.download_and_extract(force_download=force_download)
        self.extract_brca_exons(force_preparing=force_preparing)
        self.extract_sequences_bedtools(force_preparing=force_preparing)
