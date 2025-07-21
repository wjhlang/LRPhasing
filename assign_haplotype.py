#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from collections import defaultdict
from bisect import bisect_left, bisect_right
from typing import Dict, List, TextIO
import pysam

def smallest_larger(myList, pos):
    """
    Find the smallest element in a sorted list that is larger than a given position.
    """
    i = bisect_right(myList, pos)
    if i != len(myList):
        return myList[i]
    else:
        return np.inf

def assign_haplotype(read: pysam.AlignedSegment, hetsnp: pysam.VariantFile, MIN_BASE_QUALITY = 20) -> int:
    """
    Assign a haplotype to a read based on heterozygous SNPs.

    Args:
        read (pysam.AlignedSegment): The read to be phased.
        hetsnp (pysam.VariantFile): VCF file containing heterozygous SNPs.

    Returns:
        int: Haplotype (1 or 2), 3 for conflicting SNPs, or 0 for no assignment.
    """
    haplotype = 0
    # Create a dictionary of reference positions to query positions and qualities for fast lookup
    aligned_pairs = {
        ref_pos: (query_pos, read.query_qualities[query_pos])
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True)
    }

    # Fetch SNPs that overlap with the read's alignment
    for snp in hetsnp.fetch(read.reference_name, read.reference_start, read.reference_end):
        sample_data = snp.samples[0]
        if not sample_data.phased:
            continue

        snp_pos = snp.pos - 1  # VCF is 1-based, pysam is 0-based
        if snp_pos in aligned_pairs:
            query_pos, base_quality = aligned_pairs[snp_pos]

            if base_quality < MIN_BASE_QUALITY:
                continue

            base_at_snp = read.query_sequence[query_pos]
            genotype = sample_data['GT']
            
            # Determine the allele for each haplotype
            hap1_allele = snp.alleles[genotype[0]]
            hap2_allele = snp.alleles[genotype[1]]
            
            current_hap = 0
            if base_at_snp == hap1_allele:
                current_hap = 1
            elif base_at_snp == hap2_allele:
                current_hap = 2

            if current_hap != 0:
                if haplotype == 0:
                    haplotype = current_hap
                elif haplotype != current_hap:
                    return 3  # Conflicting SNPs found

    return haplotype

# --- Main Logic ---

def main():
    """Main function to parse arguments and run phasing."""
    parser = argparse.ArgumentParser(
        description='Process a BAM file to split reads by haplotype for long-read data.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--input_bam', required=True, help='Path to the input BAM file.')
    parser.add_argument('--vcf', required=True, help='Path to the phased, heterozygous SNP VCF file.')
    parser.add_argument('--output_prefix', required=True, help='Prefix for the output BAM file(s).')
    parser.add_argument(
        '--output_mode',
        choices=['separate', 'combined'],
        default='separate',
        help="Output format:\n"
             "'separate': Create different BAM files for each haplotype (hap1, hap2, unphased).\n"
             "'combined': Create a single BAM file with reads tagged with their haplotype."
    )
    args = parser.parse_args()

    # --- File Handling and Setup ---
    readfile = pysam.AlignmentFile(args.input_bam, mode='rb')
    hetsnp_vcf = pysam.VariantFile(args.vcf)

    # Sanity check: ensure VCF has exactly one sample
    if len(hetsnp_vcf.header.samples) != 1:
        print(f"Error: VCF file must contain exactly one sample. Found: {len(hetsnp_vcf.header.samples)}", file=sys.stderr)
        sys.exit(1)

    print("Caching SNP positions from VCF...")
    snp_pos_map = defaultdict(list)
    for snp in hetsnp_vcf.fetch():
        if snp.samples[0].phased:
            snp_pos_map[snp.chrom].append(snp.pos)
    for chrom in snp_pos_map:
        snp_pos_map[chrom].sort()
    print("SNP caching complete.")

    # --- Output File Initialization ---
    output_files = {}
    if args.output_mode == 'separate':
        output_files[1] = pysam.AlignmentFile(f"{args.output_prefix}.hap1.bam", "wb", template=readfile)
        output_files[2] = pysam.AlignmentFile(f"{args.output_prefix}.hap2.bam", "wb", template=readfile)
        # Unphased, conflicting, and reads with no SNP coverage
        output_files[0] = pysam.AlignmentFile(f"{args.output_prefix}.unphased.bam", "wb", template=readfile)
    else: # combined mode
        output_files[0] = pysam.AlignmentFile(f"{args.output_prefix}.phased.bam", "wb", template=readfile)

    # --- Main Processing Loop ---
    print("Processing reads...")
    processed_count = 0
    for read in readfile.fetch():
        processed_count += 1
        if processed_count % 10000 == 0:
            print(f"  ...processed {processed_count} reads")
            
        haplotype = 0 # Default for unmapped or no-SNP-overlap reads
        if not read.is_unmapped and read.reference_name in snp_pos_map:
            # Pre-check if any SNP could possibly overlap the read
            chrom_snps = snp_pos_map[read.reference_name]
            closest_snp = smallest_larger(chrom_snps, read.reference_start)
            if closest_snp <= read.reference_end:
                haplotype = assign_haplotype(read, hetsnp_vcf)

        # Set HP tag and write to the appropriate file
        if haplotype == 1:
            read.set_tag('HP', 1, value_type='i')
            writer = output_files.get(1) if args.output_mode == 'separate' else output_files[0]
        elif haplotype == 2:
            read.set_tag('HP', 2, value_type='i')
            writer = output_files.get(2) if args.output_mode == 'separate' else output_files[0]
        else: # Haplotype is 0 (unphased) or 3 (conflicting)
            if haplotype == 3:
                read.set_tag('HP', 3, value_type='i') # Tag as conflicting
            # No tag needed for 0, but it goes to the unphased/combined file
            writer = output_files[0]

        writer.write(read)

    # --- Cleanup ---
    print("Closing files...")
    for f in output_files.values():
        f.close()
    readfile.close()
    hetsnp_vcf.close()
    
    print("\n✅ Processing complete.")
    print("Output file(s) created:")
    for f_path in output_files.values():
        print(f"- {f_path.filename.decode()}")
        pysam.index(f_path.filename.decode())
        print(f"  - Index created: {f_path.filename.decode()}.bai")


if __name__ == "__main__":
    main()
