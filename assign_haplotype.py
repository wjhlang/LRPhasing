#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from collections import defaultdict, Counter
from bisect import bisect_left, bisect_right
from typing import Dict, List, TextIO, Optional, Tuple
import pysam
import math

def smallest_larger(sorted_list: List[int], pos: int) -> float:
    """
    Find the smallest element in a sorted list that is larger than a given position.
    """
    i = bisect_right(sorted_list, pos)
    return sorted_list[i] if i != len(sorted_list) else float('inf')
    
def phred_to_eps(Q: int, q_floor: int = 0, q_cap: Optional[int] = None) -> float:
    """
    Convert Phred Q to error probability ε = 10^(-Q/10), with optional floor/cap.
    Caps/floors are applied on Q before conversion.
    """
    if Q is None:
        return 0.1  # extremely conservative fallback
    if q_cap is not None and Q > q_cap:
        Q = q_cap
    if Q < q_floor:
        Q = q_floor
    return 10.0 ** (-Q / 10.0)

def log10_odds_from_eps(eps: float) -> float:
    """
    For a base that matches one haplotype allele (and mismatches the other),
    return log10((1 - eps) / eps).
    """
    # clamp eps for numeric stability and to avoid infinite/huge weights
    eps = min(max(eps, 1e-6), 0.49)
    return math.log10((1.0 - eps) / eps)

# -------------------------------
# Likelihood-based haplotype assignment
# -------------------------------

def assign_haplotype_llr(
    read: pysam.AlignedSegment,
    hetsnp: pysam.VariantFile,
    *,
    min_base_quality: int = 20,
    min_informative_snps: int = 1,
    min_log10_odds: float = 3.0,    # ~1e3:1 odds
    use_phase_set: bool = True,
    min_gq: int = 0,
    q_floor: int = 0,
    q_cap: Optional[int] = None,
) -> Tuple[int, int, float, Optional[int]]:
    """
    Assign a haplotype to a read using a log10-likelihood ratio.
    Evidence is optionally restricted to a single phase-set (PS) block.

    Returns:
        (hap, n_informative, abs_llr10, chosen_ps)

        hap:
            1 or 2  -> assigned haplotype
            3       -> conflicting (both sides seen but |LLR10| < threshold)
            0       -> unphased (too few informative SNPs or no clear evidence)

        n_informative: number of bases that matched either hap allele
        abs_llr10:     absolute log10 odds for the call (0 if unphased)
        chosen_ps:     PS value used (None if not available/used)
    """
    # Build mapping ref_pos -> (qpos, baseQ)
    aligned_pairs = {}
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if qpos is None or rpos is None:
            continue
        baseQ = 60 if read.query_qualities is None else read.query_qualities[qpos]
        aligned_pairs[rpos] = (qpos, baseQ)

    # Gather evidence per phase set
    # match: +1 if base==hap1 allele; -1 if base==hap2 allele; weight=log10 odds weight
    groups: Dict[Optional[int], List[Tuple[int, float]]] = defaultdict(list)

    if read.is_unmapped or read.reference_name is None:
        return (0, 0, 0.0, None)

    if read.query_sequence is None:
        return (0, 0, 0.0, None)

    for snp in hetsnp.fetch(read.reference_name, read.reference_start, read.reference_end):
        s = snp.samples[0]  # single-sample VCF enforced in main()
        if not s.phased:
            continue

        gt = s.get('GT')
        if not gt or len(gt) != 2 or gt[0] == gt[1]:
            # need heterozygous, phased
            continue

        if min_gq and s.get('GQ') is not None and s['GQ'] < min_gq:
            continue

        rpos0 = snp.pos - 1  # VCF 1-based -> 0-based
        if rpos0 not in aligned_pairs:
            continue

        qpos, baseQ = aligned_pairs[rpos0]
        if baseQ < min_base_quality:
            continue

        base = read.query_sequence[qpos]
        if base == 'N':
            continue

        # Alleles: genotype indices into snp.alleles
        try:
            hap1_allele = snp.alleles[gt[0]]
            hap2_allele = snp.alleles[gt[1]]
        except Exception:
            continue

        match = 0
        if base == hap1_allele:
            match = +1
        elif base == hap2_allele:
            match = -1
        else:
            # mismatch to both alleles -> uninformative for LLR here
            continue

        eps = phred_to_eps(baseQ, q_floor=q_floor, q_cap=q_cap)
        w10 = log10_odds_from_eps(eps)

        ps_val = s.get('PS') if use_phase_set else None
        groups[ps_val].append((match, w10))

    if not groups:
        return (0, 0, 0.0, None)

    # Choose which PS to use: Prefer the PS with the most informative SNPs
    best_ps = None
    best_evidence: List[Tuple[int, float]] = []
    best_key = (-1, -1.0)

    for ps, evs in groups.items():
        n = len(evs)
        tot_abs = sum(abs(w) for _, w in evs)
        key = (n, tot_abs)
        if key > best_key:
            best_key = key
            best_ps = ps
            best_evidence = evs

    informative = len(best_evidence)
    if informative < min_informative_snps:
        return (0, informative, 0.0, best_ps)

    # Compute LLR in log10 units
    llr10 = 0.0
    saw_h1 = False
    saw_h2 = False
    for m, w10 in best_evidence:
        if m == +1:
            llr10 += w10
            saw_h1 = True
        else:
            llr10 -= w10
            saw_h2 = True

    # Decision
    if llr10 >= min_log10_odds:
        return (1, informative, abs(llr10), best_ps)
    if llr10 <= -min_log10_odds:
        return (2, informative, abs(llr10), best_ps)

    # Ambiguous region where |LLR| < threshold
    if saw_h1 and saw_h2:
        return (3, informative, abs(llr10), best_ps)

    # Otherwise, unphased (|LLR| < threshold but only one haplotype seen)
    return (0, informative, abs(llr10), best_ps)

# --- Main Logic ---

def main():
    parser = argparse.ArgumentParser(
        description='Haplotype tag long reads using phased het SNPs with a likelihood ratio model.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--input_bam', required=True, help='Input BAM/CRAM (indexed).')
    parser.add_argument('--vcf', required=True, help='Phased heterozygous SNP VCF/BCF (tabix-indexed).')
    parser.add_argument('--output_prefix', required=True, help='Prefix for output BAM(s).')
    parser.add_argument('--output_mode', choices=['separate', 'combined'], default='separate',
                        help="'separate': write hap1, hap2, unphased BAMs; 'combined': one BAM with HP tags.")
    # LLR knobs
    parser.add_argument('--min_base_quality', type=int, default=20, help='Minimum per-base Q at SNP (default: 20).')
    parser.add_argument('--min_informative_snps', type=int, default=1, help='Minimum informative SNPs per read (default: 1).')
    parser.add_argument('--min_log10_odds', type=float, default=3.0, help='Min |log10 odds| to assign hap (default: 6).')
    parser.add_argument('--use_phase_set', action='store_true', default=True, help='Restrict evidence to a single PS block (default: True).')
    parser.add_argument('--no_phase_set', dest='use_phase_set', action='store_false', help='Ignore PS blocks (not recommended).')
    parser.add_argument('--min_gq', type=int, default=0, help='Ignore SNPs with FORMAT/GQ < min_gq if present (default: 0 = off).')
    parser.add_argument('--q_cap', type=int, default=None, help='Cap per-base Q at this value before ε conversion (e.g., 30 for ONT).')
    parser.add_argument('--q_floor', type=int, default=0, help='Floor per-base Q before ε conversion (default: 0).')

    args = parser.parse_args()

    readfile = pysam.AlignmentFile(args.input_bam, mode='rb')
    hetsnp_vcf = pysam.VariantFile(args.vcf)

    # Ensure VCF has exactly one sample
    if len(hetsnp_vcf.header.samples) != 1:
        print(f"Error: VCF file must contain exactly one sample. Found: {len(hetsnp_vcf.header.samples)}", file=sys.stderr)
        sys.exit(1)

    # Cache positions of phased het SNPs per chrom for a fast precheck
    print("Caching phased-het SNP positions...")
    snp_pos_map: Dict[str, List[int]] = defaultdict(list)
    for snp in hetsnp_vcf.fetch():
        s = snp.samples[0]
        if not s.phased:
            continue
        gt = s.get('GT')
        if not gt or len(gt) != 2 or gt[0] == gt[1]:
            continue
        snp_pos_map[snp.chrom].append(snp.pos)
    for chrom in snp_pos_map:
        snp_pos_map[chrom].sort()
    print("SNP caching complete.")

    # Outputs
    output_files = {}
    if args.output_mode == 'separate':
        output_files[1] = pysam.AlignmentFile(f"{args.output_prefix}.hap1.bam", "wb", template=readfile)
        output_files[2] = pysam.AlignmentFile(f"{args.output_prefix}.hap2.bam", "wb", template=readfile)
        output_files[0] = pysam.AlignmentFile(f"{args.output_prefix}.unphased.bam", "wb", template=readfile)
    else:
        output_files[0] = pysam.AlignmentFile(f"{args.output_prefix}.phased.bam", "wb", template=readfile)

    # Stats
    counts = Counter()  # keys: 'hap1','hap2','conflict','unphased','no_overlap'
    np_hist = Counter() # histogram of informative SNP count
    pc_hist = Counter() # rough histogram of confidence (abs llr10 floored)

    print("Processing reads...")
    processed = 0
    for read in readfile.fetch():
        processed += 1
        if processed % 10000 == 0:
            print(f"  ...processed {processed} reads", file=sys.stderr)

        writer = output_files[0]  # default target
        if read.is_unmapped or read.reference_name not in snp_pos_map:
            counts['no_overlap'] += 1
            read.set_tag('HP', 0, 'i')
            writer.write(read)  # goes to unphased BAM
            continue

        # Quick span check: any phased het SNP possibly overlaps?
        chrom_snps = snp_pos_map[read.reference_name]
        next_snp = smallest_larger(chrom_snps, read.reference_start)
        if next_snp > read.reference_end:
            counts['no_overlap'] += 1
            read.set_tag('HP', 0, 'i')
            writer.write(read)
            continue

        hap, n_inf, abs_llr10, ps_val = assign_haplotype_llr(
            read, hetsnp_vcf,
            min_base_quality=args.min_base_quality,
            min_informative_snps=args.min_informative_snps,
            min_log10_odds=args.min_log10_odds,
            use_phase_set=args.use_phase_set,
            min_gq=args.min_gq,
            q_floor=args.q_floor,
            q_cap=args.q_cap,
        )

        # Tag & route
        if hap == 1:
            read.set_tag('HP', 1, value_type='i')
            if ps_val is not None:
                read.set_tag('PS', int(ps_val), value_type='i')
            read.set_tag('NP', int(n_inf), value_type='i')         # informative SNP count
            read.set_tag('PC', float(abs_llr10), value_type='f')   # confidence (abs log10 odds)
            counts['hap1'] += 1
            writer = output_files.get(1) if args.output_mode == 'separate' else output_files[0]

        elif hap == 2:
            read.set_tag('HP', 2, value_type='i')
            if ps_val is not None:
                read.set_tag('PS', int(ps_val), value_type='i')
            read.set_tag('NP', int(n_inf), value_type='i')
            read.set_tag('PC', float(abs_llr10), value_type='f')
            counts['hap2'] += 1
            writer = output_files.get(2) if args.output_mode == 'separate' else output_files[0]

        elif hap == 3:
            # conflicting
            read.set_tag('HP', 3, value_type='i')
            if ps_val is not None:
                read.set_tag('PS', int(ps_val), value_type='i')
            read.set_tag('NP', int(n_inf), value_type='i')
            read.set_tag('PC', float(abs_llr10), value_type='f')
            counts['conflict'] += 1
            writer = output_files[0]  # unphased/combined file

        else:
            # unphased (no decision)
            read.set_tag('HP', 0, 'i')
            if ps_val is not None:
                read.set_tag('PS', int(ps_val), value_type='i')
            if n_inf > 0:
                read.set_tag('NP', int(n_inf), value_type='i')
                read.set_tag('PC', float(abs_llr10), value_type='f')
            counts['unphased'] += 1
            writer = output_files[0]

        np_hist[n_inf] += 1
        pc_hist[int(abs_llr10)] += 1  # coarse bin

        writer.write(read)

    # Close outputs
    print("Closing files...")
    for f in output_files.values():
        f.close()
        try:
            fname = f.filename.decode() if isinstance(f.filename, (bytes, bytearray)) else f.filename
            pysam.index(fname)
            print(f"- {fname}\n  - Index created: {fname}.bai")
        except Exception as e:
            print(f"- {f.filename}: indexing failed ({e})", file=sys.stderr)

    readfile.close()
    hetsnp_vcf.close()

    # Summary
    total = sum(counts.values())
    print("\n✅ Processing complete.")
    print("Read counts:")
    for k in ['hap1','hap2','conflict','unphased','ap']:
        v = counts[k]
        frac = (v/total*100.0) if total else 0.0
        print(f"  {k:10s}: {v:10d} ({frac:5.1f}%)")
    # Optional quick QC histos
    if np_hist:
        top_np = sorted(np_hist.items())[:10]
        print("  NP histogram (first 10 bins):", dict(top_np))
    if pc_hist:
        top_pc = sorted(pc_hist.items())[:10]
        print("  PC (|log10 odds|) histogram (first 10 bins):", dict(top_pc))

if __name__ == "__main__":
    main()
