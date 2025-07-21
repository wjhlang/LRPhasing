# LRPhasing

A simple and efficient Python script for phasing long reads from a BAM file using a VCF file of phased heterozygous SNPs. The script assigns a haplotype (HP) tag to each read based on the alleles it carries at known SNP locations.

## Features

-   Phases reads against a VCF file containing phased variants.
-   Stringent haplotype assigning, handles unphased reads, conflicting reads (spanning SNPs from both haplotypes), and reads with no SNP coverage.
-   Extremely lightweight with minimal dependencies (`pysam`).
-   Two output modes:
    1.  **`separate`**: Creates separate BAM files for Haplotype 1, Haplotype 2, and unphased/conflicting reads.
    2.  **`combined`**: Creates a single BAM file where each read is tagged with its haplotype (`HP:1`, `HP:2`, or `HP:3` for conflicting).

## Installation

1.  Clone the repository:
    ```bash
    git clone [https://github.com/YOUR_USERNAME/haplotype-phaser.git](https://github.com/YOUR_USERNAME/haplotype-phaser.git)
    cd haplotype-phaser
    ```

2.  It is recommended to use a virtual environment:
    ```bash
    python3 -m venv env
    source env/bin/activate
    ```

3.  Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

The script requires an input BAM file, a VCF file with phased SNPs, and an output prefix.

### Example Commands

**1. Separate Output Files (Default)**

This will create `results.hap1.bam`, `results.hap2.bam`, and `results.unphased.bam`.

```bash
python phase_reads.py \
    --input_bam /path/to/your/reads.bam \
    --vcf /path/to/your/phased_snps.vcf.gz \
    --output_prefix results
```

**2. Combined Output File**

This will create a single `results.phased.bam` file with `HP` tags.

```bash
python phase_reads.py \
    --input_bam /path/to/your/reads.bam \
    --vcf /path/to/your/phased_snps.vcf.gz \
    --output_prefix results \
    --output_mode combined
```

### Arguments

-   `--input_bam`: (Required) Path to the input BAM file.
-   `--vcf`: (Required) Path to the phased, heterozygous SNP VCF file (can be gzipped).
-   `--output_prefix`: (Required) Prefix for the output BAM file(s).
-   `--output_mode`: (Optional) Choose between `separate` or `combined`. Default is `separate`.

## Output

The script automatically creates BAM index (`.bai`) files for all outputs, allowing them to be loaded directly into viewers like IGV.

The `HP` tag values are:
-   `HP:1`: Read assigned to Haplotype 1.
-   `HP:2`: Read assigned to Haplotype 2.
-   `HP:3`: Read contains conflicting SNPs from both haplotypes.
-   Reads without an `HP` tag in the `unphased.bam` or `combined.bam` files had no overlapping phased SNPs.
