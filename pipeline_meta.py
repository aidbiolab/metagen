import os
import matplotlib.pyplot as plt
import pandas as pd
"""
import glob

import numpy as np
from skbio import OrdinationResults
from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

import seaborn as sns
"""
def qc_preprocessing():
    """
    Stage 1: Quality control and trimming with fastp + FastQC + MultiQC
    Assumes unpaired FASTQ in ./data/MetsNCBI/*.fastq.gz
    Uses 32 threads where possible
    """
    # Create output directories
    os.system("mkdir -p ./outputs/fastp")
    os.system("mkdir -p ./outputs/fastqc")
    os.system("mkdir -p ./outputs/multiqc")

    # Path to tools (adjust if using conda → just use binary name)
    FASTP   = "fastp"          # or just "fastp" if in PATH / conda active
    FASTQC  = "fastqc"  # or "fastqc"
    MULTIQC = "multiqc"                # or full path if using venv: "./tools/venv-multiqc/bin/multiqc"
    PIGZ    = "pigz -p 32"             # your style

    # Find all input fastq.gz (unpaired)
    input_dir = "./data/MetsNCBI"
    fastq_files = sorted(glob.glob(f"{input_dir}/*.fastq.gz") + glob.glob(f"{input_dir}/*.fq.gz"))

    if not fastq_files:
        print("No .fastq.gz or .fq.gz files found in ./data/MetsNCBI !")
        return

    print(f"Found {len(fastq_files)} samples to process.")

    # ── Run fastp + FastQC per sample ────────────────────────────────────────
    for fq in fastq_files:
        sample = os.path.basename(fq).replace(".fastq.gz", "").replace(".fq.gz", "")
        print(f"\nProcessing sample: {sample}")

        out_trim      = f"./outputs/fastp/{sample}_trimmed.fastq.gz"
        out_json      = f"./outputs/fastp/{sample}_fastp.json"
        out_html      = f"./outputs/fastp/{sample}_fastp.html"
        out_fastqc    = f"./outputs/fastqc/{sample}_fastqc"

        # fastp command (aggressive metagenomic defaults – adjust if needed)
        cmd_fastp = f"""
        {FASTP}
        -i {fq}
        -o {out_trim}
        --thread 32
        --detect_adapter_for_pe   # even for SE it tries to guess
        --cut_tail --cut_window_size 4 --cut_mean_quality 30
        --qualified_quality_phred 30 --unqualified_percent_limit 40
        --length_required 50
        -j {out_json}
        -h {out_html}
        """
        cmd_fastp = " ".join(cmd_fastp.split())
        os.system(cmd_fastp)

        # FastQC on trimmed file
        cmd_fastqc = f"""
        {FASTQC}
        {out_trim}
        -o ./outputs/fastqc/
        -t 32
        """
        cmd_fastqc = " ".join(cmd_fastqc.split())
        os.system(cmd_fastqc)

        # Optional: compress FastQC raw data if you want (html/zips are already small)
        # os.system(f"{PIGZ} ./outputs/fastqc/{sample}*fastqc.zip")

    # ── Aggregate with MultiQC ───────────────────────────────────────────────
    print("\nRunning MultiQC aggregation...")
    cmd_multiqc = f"""
    {MULTIQC}
    ./outputs/fastp ./outputs/fastqc
    -o ./outputs/multiqc
    -n metagen_qc_report
    --title "Metagenomics QC - MetsNCBI"
    --force
    """
    cmd_multiqc = " ".join(cmd_multiqc.split())
    os.system(cmd_multiqc)

    print("\nDone! Check:")
    print("  - Trimmed files:     ./outputs/fastp/*_trimmed.fastq.gz")
    print("  - Individual QC:     ./outputs/fastqc/*_fastqc.html")
    print("  - Summary report:    ./outputs/multiqc/metagen_qc_report.html")


def kraken2_classification():
    """
    Stage: Taxonomic classification with Kraken2
    Uses pre-built db: ./data/kraken/k2_standard_20260226
    Input: trimmed unpaired FASTQ from ./outputs/fastp/*_trimmed.fastq.gz
    """
    # Create output dir
    os.system("mkdir -p ./outputs/kraken2")

    # Paths
    DB_PATH = "./data/kraken/k2_standard_20260226"
    KRAKEN2_BIN = "kraken2"   # if conda; else "./tools/kraken2_bin/kraken2"

    # Check if db exists
    if not os.path.exists(f"{DB_PATH}/hash.k2d"):
        print(f"Database not found or incomplete: {DB_PATH}")
        print("Make sure you extracted the tar.gz correctly.")
        return

    # Find all trimmed fastq.gz
    input_dir = "./outputs/fastp"
    fastq_files = sorted(glob.glob(f"{input_dir}/*_trimmed.fastq.gz"))

    if not fastq_files:
        print("No trimmed FASTQ files found in ./outputs/fastp !")
        return

    print(f"Found {len(fastq_files)} trimmed samples to classify with Kraken2.")

    for fq in fastq_files:
        sample = os.path.basename(fq).replace("_trimmed.fastq.gz", "")
        print(f"\nClassifying sample: {sample}")

        out_report   = f"./outputs/kraken2/{sample}.kreport"
        out_kraken   = f"./outputs/kraken2/{sample}.kraken"
        out_classified = f"./outputs/kraken2/{sample}.classified#.fastq.gz"   # # → replaced by 1 or nothing for SE
        out_unclassified = f"./outputs/kraken2/{sample}.unclassified#.fastq.gz"

        cmd = f"""
        {KRAKEN2_BIN}
        --db {DB_PATH}
        --threads 28
        --confidence 0.1
        --gzip-compressed
        --fastq-input
        --report {out_report}
        --output {out_kraken}
        --use-names
        --report-minimizer-data
        --classified-out {out_classified}
        --unclassified-out {out_unclassified}
        {fq}
        """
        cmd = " ".join(cmd.split())   # your style: clean spaces
        print("Running:", cmd)        # optional: see command
        os.system(cmd)

    print("\nDone! Check:")
    print("  - Per-sample reports: ./outputs/kraken2/*.kreport")
    print("  - Detailed outputs:   ./outputs/kraken2/*.kraken")
    print("  - Classified reads:   ./outputs/kraken2/*.classified*.fastq.gz (if any)")
    print("\nNext: You can run Bracken on the .kreport files for abundance estimation.")

def bracken_abundance_estimation():
    read_length=50
    tax_level="S"
    threshold=10
    """
    Stage: Run Bracken on all Kraken2 .kreport files for abundance estimation.
    
    Parameters:
    - read_length: int, expected read length after trimming (e.g. 150)
    - tax_level: str, taxonomic level for abundances ('S' = species, 'G' = genus, 'F' = family, etc.)
    - threshold: int, minimum number of reads to consider a taxon (default 10)
    
    Outputs: ./outputs/bracken/*.bracken (tabular abundances per sample)
             ./outputs/bracken/*.breport (optional detailed report)
    """
    # Create output dir
    os.system("mkdir -p ./outputs/bracken")

    # Paths
    DB_PATH = "./data/kraken/k2_standard_20260226"
    REPORT_DIR = "./outputs/kraken2"

    # Check db and bracken
    if not os.path.exists(f"{DB_PATH}/hash.k2d"):
        print("Kraken DB not found!")
        return
    if not os.popen("which bracken").read().strip():
        print("Bracken not found! Install with: mamba install -c bioconda bracken -y")
        return

    # Find all .kreport files
    kreports = sorted(glob.glob(f"{REPORT_DIR}/*.kreport"))

    if not kreports:
        print("No .kreport files found in ./outputs/kraken2 !")
        return

    print(f"Found {len(kreports)} Kraken2 reports. Running Bracken...")

    for kreport in kreports:
        sample = os.path.basename(kreport).replace(".kreport", "")
        print(f"\nEstimating abundances for: {sample}")

        out_bracken = f"./outputs/bracken/{sample}.bracken"
        out_breport = f"./outputs/bracken/{sample}.breport"   # optional detailed

        cmd = f"""
        bracken
        -d {DB_PATH}
        -i {kreport}
        -o {out_bracken}
        -w {out_breport}
        -r {read_length}
        -l {tax_level}
        -t {threshold}
        """
        cmd = " ".join(cmd.split())
        print("Running:", cmd)
        os.system(cmd)

    print("\nDone! Check:")
    print(f"  - Abundance tables: ./outputs/bracken/*.bracken  (main output)")
    print(f"    Format: name | taxonomy_id | taxonomy_lvl | kraken_assigned_reads | added_reads | new_est_reads | fraction_total_reads")
    print("  - Optional reports: ./outputs/bracken/*.breport")
    print("\nNext ideas:")
    print("  - Combine all .bracken files into one table (pandas)")
    print("  - Visualize with Krona, Pavian, or Python (matplotlib/seaborn)")
    print("  - Compute alpha diversity (e.g., via KrakenTools or custom script)")
    print("  - Filter human/contaminant reads if needed")

def merge_bracken_outputs():
    value_column="fraction_total_reads"
    output_file="./outputs/bracken/merged_rel_abund.csv"
    #merge_bracken_outputs(value_column="new_est_reads", output_file="./outputs/bracken/merged_counts.csv")
    """
    Merge all .bracken files into one table (robust against column overlap).
    """
    bracken_dir = "./outputs/bracken"
    bracken_files = sorted(glob.glob(f"{bracken_dir}/*.bracken"))

    if not bracken_files:
        print("No .bracken files found!")
        return

    print(f"Found {len(bracken_files)} Bracken files. Merging...")

    # ── Step 1: Collect only abundance columns ───────────────────────────────
    abundance_dfs = []
    for file in bracken_files:
        sample = os.path.basename(file).replace(".bracken", "")
        print(f"Processing: {sample}")

        df = pd.read_csv(file, sep="\t", usecols=["name", value_column])
        df = df.rename(columns={value_column: sample})
        df.set_index("name", inplace=True)
        abundance_dfs.append(df)

    # Outer merge on taxon name (index)
    merged = pd.concat(abundance_dfs, axis=1, join="outer").fillna(0)

    # ── Step 2: Add taxonomy metadata from FIRST file only ───────────────────
    first_file = bracken_files[0]
    tax_df = pd.read_csv(first_file, sep="\t", usecols=["name", "taxonomy_id", "taxonomy_lvl"])
    tax_df.set_index("name", inplace=True)

    # Join with suffixes to handle any potential overlap safely
    merged = merged.join(tax_df, how="left", rsuffix="_tax")  # rsuffix prevents overlap error

    # If suffix was added (rare), clean up column names
    if any(col.endswith("_tax") for col in merged.columns):
        print("Warning: taxonomy columns had overlap → using _tax suffix")
        merged = merged.rename(columns=lambda c: c.replace("_tax", "") if c.endswith("_tax") else c)

    # Reorder: taxonomy first, then samples
    tax_cols = ["taxonomy_id", "taxonomy_lvl"]
    sample_cols = [col for col in merged.columns if col not in tax_cols]
    merged = merged[tax_cols + sample_cols]

    # Optional: sort by total abundance
    merged["total_est_reads"] = merged[sample_cols].sum(axis=1)
    merged = merged.sort_values("total_est_reads", ascending=False).drop(columns="total_est_reads")

    # Save
    merged.to_csv(output_file)
    print(f"\nMerged table saved to: {output_file}")
    print(f"Shape: {merged.shape} (taxa × samples)")
    print("Top 5 taxa preview:")
    print(merged.head().to_string())

def generate_krona_per_sample():
    """
    Generate interactive Krona HTML charts for each Kraken2 .kreport file.
    Outputs go to ./outputs/krona/
    """
    os.system("mkdir -p ./outputs/krona")

    kreport_files = sorted(glob.glob("./outputs/kraken2/*.kreport"))

    if not kreport_files:
        print("No .kreport files found in ./outputs/kraken2/ !")
        print("Make sure Kraken2 finished successfully.")
        return

    print(f"Found {len(kreport_files)} Kraken2 reports. Generating Krona charts...")

    for kreport in kreport_files:
        sample = os.path.basename(kreport).replace(".kreport", "")
        out_html = f"./outputs/krona/{sample}.krona.html"

        # Correct command for Kraken2 .kreport files:
        # -t 5 : taxonomy ID is in column 5 (standard in .kreport)
        # -m 3 : magnitude (reads assigned) is in column 3
        cmd = f"""
        ktImportTaxonomy
        -t 5
        -m 3
        -o {out_html}
        {kreport}
        """

        cmd = " ".join(cmd.split())
        print(f"Running for {sample}: {cmd}")
        os.system(cmd)

    print("\nDone! Krona charts created in:")
    print("  ./outputs/krona/*.krona.html")
    print("\nTo view:")
    print("  - Download the .html files to your local computer (via scp, FileZilla, etc.)")
    print("  - Open in any modern browser (Chrome/Firefox/Edge)")
    print("  - Interactive sunburst chart — click sectors to zoom into taxa")


def compute_alpha_diversity():

    input_csv="./outputs/bracken/merged_rel_abund.csv"
    metrics=[
        "shannon",
        "simpson",
        "enspie",
        "pielou_evenness",
        "observed_features",
        "chao1",
        "berger_parker_d"]
    output_csv="./outputs/diversity/alpha_diversity_counts.csv"

    df = pd.read_csv(input_csv, index_col=0)

    abund_cols = [c for c in df.columns if c not in ["taxonomy_id", "taxonomy_lvl"]]
    abund = df[abund_cols].T  # samples as rows, taxa as columns
    data = abund.values.astype(float)
    sample_ids = abund.index.tolist()
    
    results = {}
    for metric in metrics:
        try:
            results[metric] = alpha_diversity(metric, data, ids=sample_ids)
            print(f"Computed: {metric}")
        except ValueError as e:
            print(f"Skipped '{metric}': {e}")
    
    alpha_df = pd.DataFrame(results)
    alpha_df.index.name = "sample"
    
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    alpha_df.to_csv(output_csv)
    
    print(f"\nAlpha diversity saved to: {output_csv}")
    print(f"Shape: {alpha_df.shape}")
    print("Preview (first 10 samples):")
    print(alpha_df.head(10).to_string())
    
    return alpha_df

def compute_beta_diversity():
    input_csv="./outputs/bracken/merged_rel_abund.csv"
    metric="braycurtis"
    output_pcoa_csv="./outputs/diversity/beta_pcoa_braycurtis_counts.csv"

    """
    Compute beta diversity distance matrix and PCoA ordination coordinates.
    Saves only the PCoA coordinates (no plot generated).
    """
    
    df = pd.read_csv(input_csv, index_col=0)
    
    # Select abundance columns only
    abund_cols = [c for c in df.columns if c not in ["taxonomy_id", "taxonomy_lvl"]]
    abund = df[abund_cols].T  # samples × taxa
    
    data = abund.values.astype(float)
    samples = abund.index.tolist()
    
    print(f"Computing {metric} distance matrix for {len(samples)} samples...")
    
    dm = beta_diversity(metric, data, ids=samples)
    
    print("Running PCoA ordination...")
    ordination = pcoa(dm)
    
    # Extract first 3 PCs
    coords = ordination.samples.iloc[:, :3]
    coords.columns = ["PC1", "PC2", "PC3"]
    coords.index.name = "sample"
    
    # Explained variance – use .iloc for positional access
    explained = ordination.proportion_explained
    print("\nExplained variance:")
    print(f"PC1: {explained.iloc[0]:.3%}")
    print(f"PC2: {explained.iloc[1]:.3%}")
    print(f"PC3: {explained.iloc[2]:.3%}")
    
    # Save coordinates
    os.makedirs("./outputs/diversity", exist_ok=True)
    coords.to_csv(output_pcoa_csv)
    
    print(f"\nSaved to: {output_pcoa_csv}")
    print("Shape:", coords.shape)
    print("Preview (first 5 samples):")
    print(coords.head().to_string())
    
    return coords


def plot_stacked_bar_per_sample_top(
    input_csv="./outputs/bracken/merged_rel_abund.csv",
    top_n_per_sample=6,
    output_file="./outputs/diversity/stacked_bar_per_sample_top.png",
    figsize=(16, 9)
):
    df = pd.read_csv(input_csv, index_col=0)
    sample_cols = [c for c in df.columns if c not in ["taxonomy_id", "taxonomy_lvl"]]
    abund = df[sample_cols].copy() * 100  # to %

    fig, ax = plt.subplots(figsize=figsize)
    bottom = pd.Series(0, index=sample_cols)

    # Sort samples by some criteria (optional: e.g. by Shannon or dominant taxon)
    sample_order = sample_cols  # or sort by abund.sum() or alpha diversity

    for i in range(top_n_per_sample):
        # Get current top taxon across remaining abundance
        remaining = abund - bottom
        current_top = remaining.sum(axis=0).idxmax()
        if remaining[current_top].sum() == 0:
            break

        ax.bar(
            sample_order,
            remaining[current_top][sample_order],
            bottom=bottom[sample_order],
            label=current_top,
            width=0.8
        )
        bottom += remaining[current_top]

    # Add "Other" for the rest
    other = 100 - bottom
    ax.bar(sample_order, other[sample_order], bottom=bottom[sample_order], label="Other", color="lightgray")

    ax.set_ylim(0, 100)
    ax.set_ylabel("Relative Abundance (%)")
    ax.set_title(f"Top {top_n_per_sample} Taxa per Sample + Other")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", title="Taxa")
    plt.xticks(rotation=90, ha="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Per-sample top taxa plot saved to: {output_file}")

#qc_preprocessing()
#kraken2_classification()
#bracken_abundance_estimation()
#generate_krona_per_sample()
#merge_bracken_outputs()
#compute_alpha_diversity()
#compute_beta_diversity()
plot_stacked_bar_per_sample_top(top_n_per_sample=5)
