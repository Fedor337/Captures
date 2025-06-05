import argparse
from reference_preparer import ReferencePreparer
from probe_generator import ProbeGenerator

def main():
    """
    Command-line pipeline for designing BRCA1/2 hybridization probes.

    This script performs the following steps:
    1. Downloads and prepares genome and annotation data
    2. Extracts exon coordinates for BRCA1/2
    3. Generates overlapping probes of specified length and step

    Usage:
        python main.py [--force-download] [--force-prep] \
                       [--input-fasta path] [--output-fasta path] \
                       [--probe-length N] [--max-step N]

    Example:
        python main.py --force-download --probe-length 120 --max-step 50
    """
    parser = argparse.ArgumentParser(description="BRCA1/2 Probe Design Pipeline")
    parser.add_argument("--force-download", action="store_true",
                        help="Force re-download and re-extraction of reference files")
    parser.add_argument("--force-prep", action="store_true",
                        help="Force regeneration of BED/FASTA exon data")
    parser.add_argument("--input-fasta", type=str, default="data/brca_exons.fa",
                        help="Path to exon FASTA file (default: data/brca_exons.fa)")
    parser.add_argument("--output-fasta", type=str, default="data/brca_probes.fa",
                        help="Path to output probe FASTA file (default: data/brca_probes.fa)")
    parser.add_argument("--probe-length", type=int, default=120,
                        help="Length of each probe (default: 120)")
    parser.add_argument("--max-step", type=int, default=60,
                        help="Maximum step size between probes (default: 60)")

    args = parser.parse_args()

    print("[1] Preparing reference data...")
    rp = ReferencePreparer()
    rp.prepare_all(force_download=args.force_download, force_preparing=args.force_prep)

    print("[2] Generating probes...")
    pg = ProbeGenerator(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        probe_length=args.probe_length,
        max_step=args.max_step
    )
    pg.generate_all()

    print("[OK] Pipeline complete.")

if __name__ == "__main__":
    main()
