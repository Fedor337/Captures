import argparse
from reference_preparer import ReferencePreparer
from probe_generator import ProbeGenerator
from probe_filter_pipeline import ProbeFilterPipeline

def main():
    """
    Command-line pipeline for designing BRCA1/2 hybridization probes.

    This script performs the following steps:
    1. Downloads and prepares genome and annotation data
    2. Extracts exon coordinates for BRCA1/2
    3. Generates overlapping probes of specified length and step
    4. Applies optional filters: GC content, Tm, repeats, structure

    Usage:
        python main.py [--force-download] [--force-prep] \
                       [--input-fasta path] [--output-fasta path] \
                       [--probe-length N] [--max-step N] \
                       [--gc-min N] [--gc-max N] \
                       [--tm-min N] [--tm-max N] \
                       [--no-repeats] [--structure-filter] [--dg-threshold N]

    Example:
        python main.py --force-download --probe-length 120 --max-step 60 \
                       --gc-min 40 --gc-max 60 --tm-min 65 --tm-max 72 \
                       --no-repeats --structure-filter --dg-threshold -9.0
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

    # Filter options
    parser.add_argument("--gc-min", type=float, default=40.0, help="Minimum GC content (%)")
    parser.add_argument("--gc-max", type=float, default=60.0, help="Maximum GC content (%)")
    parser.add_argument("--tm-min", type=float, default=65.0, help="Minimum melting temperature (°C)")
    parser.add_argument("--tm-max", type=float, default=72.0, help="Maximum melting temperature (°C)")
    parser.add_argument("--no-repeats", action="store_true", help="Remove probes with repeats")
    parser.add_argument("--structure-filter", action="store_true", help="Enable secondary structure filtering via RNAfold")
    parser.add_argument("--dg-threshold", type=float, default=-9.0, help="Minimum acceptable ΔG (kcal/mol) for RNAfold")

    args = parser.parse_args()

    print("[1] Preparing reference data...")
    rp = ReferencePreparer()
    rp.prepare_all(force_download=args.force_download, force_preparing=args.force_prep)

    print("[2] Generating probes...")
    intermediate_fasta = args.output_fasta.replace(".fa", ".raw.fa")
    pg = ProbeGenerator(
        input_fasta=args.input_fasta,
        output_fasta=intermediate_fasta,
        probe_length=args.probe_length,
        max_step=args.max_step
    )
    pg.generate_all()

    print("[3] Applying filters...")
    pf = ProbeFilterPipeline(
        input_fasta=intermediate_fasta,
        output_fasta=args.output_fasta,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        tm_min=args.tm_min,
        tm_max=args.tm_max,
        allow_repeats=not args.no_repeats,
        structure_filter=args.structure_filter,
        dg_threshold=args.dg_threshold
    )
    pf.apply_all()

    print("[OK] Pipeline complete.")

if __name__ == "__main__":
    main()
