# LRPhasing

Lightweight long-read phasing with a **log-likelihood-ratio (LLR)** model.  
Given an indexed BAM/CRAM and a tabix-indexed VCF of *phased, heterozygous* SNPs (single sample), LRPhasing assigns each read to a haplotype and writes either per-haplotype BAMs or a single BAM tagged with haplotype and confidence.

- **Minimal deps:** `pysam` (htslib under the hood)  
- **Robust LLR:** per-base Phred Q → error ε → weight `log10((1−ε)/ε)`; evidence optionally restricted to a single phase-set (`PS`)  
- **Clear tagging:** `HP` (haplotype), `NP` (informative SNP count), `PC` (|log10-odds|), and optional `PS` (phase-set used)

---

## What’s new (LLR model)

- Replaced simple counting with a **log10-odds** decision rule.  
- Supports **quality capping/flooring** (`--q_cap`, `--q_floor`) for platform-aware weighting (e.g., ONT).  

LLR per SNP uses base error probability from Phred Q:  
`ε = 10^(-Q/10)`, weight `w = log10((1 − ε) / ε)`.

---

## Installation

```bash
git clone https://github.com/wjhlang/LRPhasing.git
cd LRPhasing

# (recommended) use a virtualenv
python3 -m venv env
source env/bin/activate

pip install -r requirements.txt
```

**Input requirements**

- BAM must be indexed (`.bai`)  
- VCF/BCF must be **single-sample**, **phased, heterozygous SNPs**, bgzipped and **tabix-indexed** (`.vcf.gz` + `.tbi`)

---

## Usage

```bash
python phase_reads.py \
  --input_bam /path/to/reads.bam \
  --vcf /path/to/phased_snps.vcf.gz \
  --output_prefix results
```

### Output modes

**1) Separate (default)**  
Writes:
- `results.hap1.bam`
- `results.hap2.bam`
- `results.unphased.bam` (contains **unphased + conflicting + no-overlap** reads)

```bash
python phase_reads.py \
  --input_bam /path/to/reads.bam \
  --vcf /path/to/phased_snps.vcf.gz \
  --output_prefix results
```

**2) Combined**  
Writes a single `results.phased.bam` with per-read tags.

```bash
python phase_reads.py \
  --input_bam /path/to/reads.bam \
  --vcf /path/to/phased_snps.vcf.gz \
  --output_prefix results \
  --output_mode combined
```

All outputs are indexed automatically.

---

## Arguments (current defaults)

**Core**
- `--input_bam` *(required)*: BAM/CRAM (indexed)  
- `--vcf` *(required)*: phased, het SNP VCF/BCF (tabix-indexed)  
- `--output_prefix` *(required)*  
- `--output_mode {separate,combined}`: default **separate**

**LLR knobs**
- `--min_base_quality INT` (default **20**): per-base Q threshold at SNPs  
- `--min_informative_snps INT` (default **1**): minimum matching SNPs used to attempt a call  
- `--min_log10_odds FLOAT` (default **3.0**): require \|LLR10\| ≥ threshold to assign HP=1/2  
  *(3.0 ≈ 10³:1 odds; raise for stricter calls)*  
- `--use_phase_set / --no_phase_set` (default **use_phase_set**): restrict evidence to a single `PS` block  
- `--min_gq INT` (default **0**): if present in FORMAT, ignore SNPs with GQ < value  
- `--q_cap INT` (default **None**): cap Q before ε conversion (e.g., `--q_cap 30` for ONT)  
- `--q_floor INT` (default **0**): floor Q before ε conversion

> The script **enforces exactly one sample** in the VCF and will exit otherwise.

---

## Read tags & files

In **combined** mode, every read goes to `results.phased.bam`.  
In **separate** mode, routing depends on outcome.

| Case | Tags |
|---|---|
| **Haplotype 1** | `HP:i:1`, `NP:i:n`, `PC:f:x`, optional `PS:i:val` |
| **Haplotype 2** | `HP:i:2`, `NP:i:n`, `PC:f:x`, optional `PS:i:val` |
| **Conflicting** | `HP:i:3`, `NP:i:n`, `PC:f:x`, optional `PS:i:val` |
| **Unphased (weak/insufficient evidence)** | **`HP:i:0`**, optional `PS:i:val`; `NP`/`PC` may be present if any evidence was considered |
| **No-overlap (no SNP overlapped at all)** | **`HP:i:0` only** (no `NP`/`PC`/`PS`) |

**Distinguishing `unphased` vs `no-overlap` when both have `HP:0`:**
- Likely **no-overlap**: `HP:0` and **no other tags** (`NP`/`PC`/`PS` absent).
- Likely **unphased**: `HP:0` **and** `NP` and/or `PC` present (and sometimes `PS`).

**Summary** printed to STDOUT includes counts for `hap1`, `hap2`, `conflict`, `unphased`, and `no_overlap`, plus quick histograms for `NP` and `PC`.

---

## Examples

**ONT with conservative quality weighting**
```bash
python phase_reads.py \
  --input_bam ont.bam \
  --vcf sample.phased.vcf.gz \
  --output_prefix ont_llr \
  --q_cap 30 \
  --min_base_quality 15 \
  --min_log10_odds 4.0
```

## Notes & best practices

- **VCF must be phased & het**: unphased or homozygous records are ignored.  
- **Single-sample VCF only**: script exits if multiple samples are present.  
- **Conflicting vs unphased**:  
  - *Conflicting* (`HP=3`): both hap1 and hap2 evidence present but below threshold.  
  - *Unphased*: Insufficient or one-sided weak evidence (\|LLR\| < threshold).  
- **Platform tuning**: `--q_cap` helps prevent overweighting inflated per-base Q (common on ONT).  
- **Thresholds**: Increase `--min_log10_odds` and/or `--min_informative_snps` to trade recall for precision.  
- **CRAM**: Supported via pysam; ensure reference is accessible (`REF_PATH`, `.fai`, or environment config).

---

## License

MIT
