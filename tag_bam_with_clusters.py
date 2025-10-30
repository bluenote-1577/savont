#!/usr/bin/env python3
"""
Tag BAM file reads with HP:i:NUM based on cluster assignments from final_clusters.tsv

Usage:
    python tag_bam_with_clusters.py input.bam final_clusters.tsv output.bam
"""

import sys
import argparse
import pysam
from pathlib import Path


def parse_cluster_file(cluster_file):
    """
    Parse the final_clusters.tsv file and return a mapping of read_id -> cluster_number

    Args:
        cluster_file: Path to final_clusters.tsv

    Returns:
        dict: Mapping of read_id (UUID) to cluster number
    """
    read_to_cluster = {}
    current_cluster_num = None

    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Check if this is a cluster header line
            if line.startswith('final_cluster_'):
                # Extract cluster number from header like "final_cluster_573"
                parts = line.split('\t')
                cluster_name = parts[0]
                current_cluster_num = int(cluster_name.split('_')[-1])
            else:
                # This is a member line with format: "uuid similarity_score"
                parts = line.split()
                if len(parts) >= 1:
                    read_id = parts[0]
                    if current_cluster_num is not None:
                        read_to_cluster[read_id] = current_cluster_num

    return read_to_cluster


def tag_bam_file(input_bam, cluster_mapping, output_bam):
    """
    Read input BAM, add HP:i:NUM tags based on cluster mapping, write to output BAM

    Args:
        input_bam: Path to input BAM file
        cluster_mapping: dict mapping read_id to cluster number
        output_bam: Path to output BAM file
    """
    stats = {
        'total_reads': 0,
        'tagged_reads': 0,
        'untagged_reads': 0
    }

    with pysam.AlignmentFile(input_bam, 'rb') as bam_in:
        with pysam.AlignmentFile(output_bam, 'wb', header=bam_in.header) as bam_out:
            for read in bam_in:
                stats['total_reads'] += 1

                # Get the read name (should match the UUID in the cluster file)
                read_name = read.query_name

                # Check if this read is in a cluster
                if read_name in cluster_mapping:
                    cluster_num = cluster_mapping[read_name]
                    # Add HP tag as integer
                    read.set_tag('HP', cluster_num, value_type='i')
                    stats['tagged_reads'] += 1
                else:
                    stats['untagged_reads'] += 1

                # Write the read to output BAM
                bam_out.write(read)

    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Tag BAM file reads with HP:i:NUM based on cluster assignments'
    )
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('cluster_file', help='Cluster file (final_clusters.tsv)')
    parser.add_argument('output_bam', help='Output BAM file with HP tags')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print verbose output')

    args = parser.parse_args()

    # Validate input files exist
    if not Path(args.input_bam).exists():
        print(f"Error: Input BAM file not found: {args.input_bam}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.cluster_file).exists():
        print(f"Error: Cluster file not found: {args.cluster_file}", file=sys.stderr)
        sys.exit(1)

    # Parse cluster file
    if args.verbose:
        print(f"Parsing cluster file: {args.cluster_file}")

    cluster_mapping = parse_cluster_file(args.cluster_file)

    if args.verbose:
        print(f"Found {len(cluster_mapping)} reads in {len(set(cluster_mapping.values()))} clusters")

    # Tag BAM file
    if args.verbose:
        print(f"Processing BAM file: {args.input_bam}")

    stats = tag_bam_file(args.input_bam, cluster_mapping, args.output_bam)

    # Print statistics
    print(f"BAM tagging complete!")
    print(f"  Total reads processed: {stats['total_reads']}")
    print(f"  Reads tagged with HP: {stats['tagged_reads']}")
    print(f"  Reads without cluster: {stats['untagged_reads']}")
    print(f"  Output written to: {args.output_bam}")

    # Index the output BAM if possible
    if args.verbose:
        print(f"Indexing output BAM file...")

    try:
        pysam.index(args.output_bam)
        if args.verbose:
            print(f"Index created: {args.output_bam}.bai")
    except Exception as e:
        print(f"Warning: Could not index output BAM: {e}", file=sys.stderr)


if __name__ == '__main__':
    main()
