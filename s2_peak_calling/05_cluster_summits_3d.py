#!/usr/bin/env python3
"""
3D Peak Clustering with MM (Distance-Weighted) Centroid

This script clusters genomic peaks using 3D euclidean distance of
(start, summit, end) and calculates centroids using MM (Majorize-Minimization)
distance-weighted algorithm.

=== Clustering ===
• 3D distance: sqrt((start_i - start_j)^2 + (summit_i - summit_j)^2 + (end_i - end_j)^2)
• Adjacent peaks (sorted by summit) within max_gap 3D distance are grouped
• Peak shape (length, asymmetry) naturally affects clustering

=== Centroid Calculation (MM Distance-Weighted) ===
• Single peaks: Original coordinates preserved
• Clustered peaks: MM distance-weighted centroid with asymmetric boundaries
  - Summit: MM algorithm (weight = 1/distance, iterative convergence)
  - Left boundary: MM distance-weighted calculation of (summit - start) lengths
  - Right boundary: MM distance-weighted calculation of (end - summit) lengths

=== Output Files ===
1. Main centroid file (always): Merged peaks representing each cluster
2. Sample-specific files (optional with --save_clustered_peaks):
   Original peaks grouped by sample, each with cluster_id for tracking
3. Statistics files (optional):
   - total.tsv (with --save_total): Per-cluster quality metrics
   - summary.tsv (with --save_summary): Overall quality metrics

=== Output Format ===
All output files use narrowPeak format (10 columns, tab-separated, no header):
  Column 1-5: chr, start, end, name, score
  Column 6: cluster_id (e.g., 'chr1_0') - replaces strand for tracking
  Column 7-10: signalValue, pValue, qValue, peak_offset
"""
import pandas as pd
import numpy as np
import argparse
from scipy.spatial.distance import pdist
from multiprocessing import Pool
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')


def create_3d_coordinates(start_vals, summit_vals, end_vals):
    """Create 3D coordinate array from peak positions."""
    return np.column_stack([start_vals, summit_vals, end_vals])


def cluster_peaks_3d(start_vals, summit_vals, end_vals, max_gap=70):
    """
    Cluster peaks based on 3D euclidean distance of (start, summit, end).

    Peaks are sorted by summit position, then adjacent peaks within
    max_gap 3D distance are grouped together. This naturally incorporates
    peak shape (length, asymmetry) into clustering.

    Args:
        start_vals: Array of peak start positions
        summit_vals: Sorted array of summit positions
        end_vals: Array of peak end positions
        max_gap: Maximum 3D euclidean distance allowed within a cluster

    Returns:
        List of lists, each containing indices of peaks in a cluster
    """
    if len(summit_vals) == 0:
        return []

    clusters = []
    current = [0]  # Store indices

    for i in range(1, len(summit_vals)):
        dist = np.sqrt(
            (start_vals[i] - start_vals[i-1]) ** 2 +
            (summit_vals[i] - summit_vals[i-1]) ** 2 +
            (end_vals[i] - end_vals[i-1]) ** 2
        )
        if dist <= max_gap:
            current.append(i)
        else:
            clusters.append(current)
            current = [i]
    clusters.append(current)

    return clusters


def mm_centroid_1d(values, max_iter=100, tol=1e-6):
    """
    Calculate MM (Majorize-Minimization) distance-weighted centroid for 1D values.

    Uses iterative algorithm where weight = 1/distance.
    Closer points have higher influence on the centroid.
    This approximates the geometric median (L1 centroid).

    Args:
        values: Array of 1D positions
        max_iter: Maximum iterations
        tol: Convergence tolerance

    Returns:
        Converged centroid position
    """
    if len(values) == 1:
        return values[0]

    c = np.mean(values)
    for _ in range(max_iter):
        distances = np.abs(values - c) + 1e-6  # Avoid division by zero
        weights = 1.0 / distances
        new_c = np.average(values, weights=weights)
        if abs(new_c - c) < tol:
            break
        c = new_c
    return c


def create_centroid(start_vals, summit_vals, end_vals, sample_vals, chrom,
                    raw_score_vals, signal_vals, pvalue_vals, qvalue_vals, cluster_id="."):
    """
    Create centroid from cluster peaks using MM distance-weighted algorithm.

    Single peak (n=1): Preserves original coordinates
    Multiple peaks (n>1): MM distance-weighted centroid with asymmetric boundaries
        - Weighted summit: MM algorithm (weight = 1/distance, iterative)
        - Left/Right lengths: MM distance-weighted from input peaks

    Returns:
        List representing narrowPeak format row
    """
    if len(summit_vals) == 1:
        # Single peak: preserve original coordinates
        start = int(start_vals[0])
        end = int(end_vals[0])
        summit = int(summit_vals[0])
        peak_offset = summit - start
        samples = str(sample_vals[0])
        sum_raw_score = round(float(raw_score_vals[0]), 4)
        sum_signal = round(float(signal_vals[0]), 4)
        sum_pvalue = round(float(pvalue_vals[0]), 4)
        sum_qvalue = round(float(qvalue_vals[0]), 4)
        return [chrom, start, end, samples, sum_raw_score, cluster_id,
                sum_signal, sum_pvalue, sum_qvalue, peak_offset]

    # Multi-peak: calculate MM distance-weighted centroid
    # Summit: MM algorithm (weight = 1/distance)
    weighted_summit = int(round(mm_centroid_1d(summit_vals.astype(float))))

    # Calculate left and right lengths from each input peak
    left_lengths = summit_vals - start_vals
    right_lengths = end_vals - summit_vals

    # Get MM distance-weighted values for left and right lengths
    left_mm = int(round(mm_centroid_1d(left_lengths.astype(float))))
    right_mm = int(round(mm_centroid_1d(right_lengths.astype(float))))

    # Create asymmetric peak around weighted summit
    start = weighted_summit - left_mm
    if start < 0:
        start = 0
    end = weighted_summit + right_mm
    peak_offset = left_mm
    samples = ",".join(sorted(set(sample_vals)))

    # Sum all score columns
    sum_raw_score = round(float(np.sum(raw_score_vals)), 4)
    sum_signal = round(float(np.sum(signal_vals)), 4)
    sum_pvalue = round(float(np.sum(pvalue_vals)), 4)
    sum_qvalue = round(float(np.sum(qvalue_vals)), 4)

    return [chrom, start, end, samples, sum_raw_score, cluster_id,
            sum_signal, sum_pvalue, sum_qvalue, peak_offset]


def calculate_distance_metrics(distances, max_gap):
    """Calculate distance-based cluster quality metrics."""
    if len(distances) == 0:
        return {'mean_dist': 0.0, 'std_dist': 0.0, 'compactness': 0.0, 'cv': 0.0, 'ratio': 0.0}

    mean_dist = np.mean(distances)
    std_dist = np.std(distances, ddof=1) if len(distances) > 1 else 0.0
    cv = std_dist / mean_dist if mean_dist > 0 else 0.0
    ratio = mean_dist / max_gap if max_gap > 0 else 0.0

    return {
        'mean_dist': mean_dist,
        'std_dist': std_dist,
        'compactness': mean_dist,
        'cv': cv,
        'ratio': ratio
    }


def collect_distance_stats(positions_3d, length_vals, cluster_id, max_gap):
    """Collect distance-based cluster quality metrics."""
    n = len(positions_3d)

    if n < 1:
        return {
            'cluster_id': cluster_id,
            'n_peaks': n,
            'mean_dist': 0.0,
            'std_dist': 0.0,
            'compactness': 0.0,
            'cv_dist': 0.0,
            'ratio': 0.0,
            'mean_length': 0.0,
            'cv_length': 0.0
        }

    # Distance metrics (using 3D coordinates)
    if n >= 2:
        distances = pdist(positions_3d, metric='euclidean')
        dist_metrics = calculate_distance_metrics(distances, max_gap)
    else:
        dist_metrics = {'mean_dist': 0.0, 'std_dist': 0.0, 'compactness': 0.0, 'cv': 0.0, 'ratio': 0.0}

    # Length metrics
    mean_length = np.mean(length_vals)
    std_length = np.std(length_vals, ddof=1) if n > 1 else 0.0
    cv_length = std_length / mean_length if mean_length > 0 else 0.0

    return {
        'cluster_id': cluster_id,
        'n_peaks': n,
        'mean_dist': dist_metrics['mean_dist'],
        'std_dist': dist_metrics['std_dist'],
        'compactness': dist_metrics['compactness'],
        'cv_dist': dist_metrics['cv'],
        'ratio': dist_metrics['ratio'],
        'mean_length': mean_length,
        'cv_length': cv_length
    }


def process_chromosome(chrom_data):
    """Process single chromosome with 3D max_gap-based clustering."""
    chrom, df_chr, max_gap, min_cluster_size, save_total, save_summary, save_clustered_peaks = chrom_data

    if len(df_chr) == 0:
        return [], 0, 0, [], [], []

    # Sort by summit position
    df_chr = df_chr.sort_values('pos').reset_index(drop=True)

    # Cluster using 3D distance (start, summit, end)
    cluster_indices = cluster_peaks_3d(
        df_chr['start'].values, df_chr['pos'].values, df_chr['end'].values,
        max_gap=max_gap
    )

    centroid_rows = []
    total_clusters = len(cluster_indices)
    valid_clusters = 0
    distance_stats = []
    all_distances = []
    all_peaks_list = []

    # Process each cluster
    for cluster_id, indices in enumerate(cluster_indices):
        if len(indices) < min_cluster_size:
            continue

        valid_clusters += 1
        sub = df_chr.iloc[indices]

        start_vals = sub["start"].values
        summit_vals = sub["pos"].values
        end_vals = sub["end"].values
        length_vals = sub["length"].values
        sample_vals = sub["sample"].values
        raw_score_vals = sub["raw_score"].values
        signal_vals = sub["signalValue"].values
        pvalue_vals = sub["pValue"].values
        qvalue_vals = sub["qValue"].values

        cluster_label = f"{chrom}_{cluster_id}"

        # Save original peaks if requested
        if save_clustered_peaks:
            for idx in indices:
                row = df_chr.iloc[idx]
                all_peaks_list.append({
                    'chr': chrom,
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'sample': str(row['sample']),
                    'raw_score': float(row['raw_score']),
                    'cluster_id': cluster_label,
                    'signalValue': float(row['signalValue']),
                    'pValue': float(row['pValue']),
                    'qValue': float(row['qValue']),
                    'peak_offset': int(row['peak_offset'])
                })

        # Calculate 3D coordinates for statistics
        cluster_positions_3d = create_3d_coordinates(start_vals, summit_vals, end_vals)

        if len(cluster_positions_3d) >= 2:
            cluster_distances_3d = pdist(cluster_positions_3d, metric='euclidean')
            if save_summary:
                all_distances.extend(cluster_distances_3d)

        # Collect statistics
        if save_total:
            stats = collect_distance_stats(cluster_positions_3d, length_vals, cluster_label, max_gap)
            stats['chromosome'] = chrom
            distance_stats.append(stats)

        # Create centroid
        centroid_row = create_centroid(
            start_vals, summit_vals, end_vals, sample_vals, chrom,
            raw_score_vals, signal_vals, pvalue_vals, qvalue_vals, cluster_label
        )
        centroid_rows.append(centroid_row)

    return centroid_rows, total_clusters, valid_clusters, distance_stats, all_distances, all_peaks_list


def save_clustered_peaks_by_sample(all_clustered_peaks, output_file):
    """
    Save original peaks grouped by sample to separate files.

    File naming: {output_prefix}_{sample}.narrowPeak
    """
    if len(all_clustered_peaks) == 0:
        return []

    peaks_df = pd.DataFrame(all_clustered_peaks)
    output_prefix = output_file.rsplit('.', 1)[0]
    saved_files = []

    for sample in sorted(peaks_df['sample'].unique()):
        sample_peaks = peaks_df[peaks_df['sample'] == sample].copy()
        sample_peaks = sample_peaks.sort_values(['chr', 'start'])

        output_data = sample_peaks[['chr', 'start', 'end', 'sample', 'raw_score',
                                     'cluster_id', 'signalValue', 'pValue', 'qValue', 'peak_offset']]

        sample_output = f"{output_prefix}_{sample}.narrowPeak"
        output_data.to_csv(sample_output, sep='\t', header=False, index=False)
        saved_files.append((sample, sample_output, len(sample_peaks)))

    return saved_files


def save_distance_statistics(distance_stats_df, output_file, all_distances, save_total=True, save_summary=True):
    """Save distance-based cluster quality statistics to TSV files."""
    saved_files = []

    # Save per-cluster statistics
    if save_total and distance_stats_df is not None:
        distance_stats_df = distance_stats_df.fillna(0)
        total_file = output_file.rsplit('.', 1)[0] + '_total.tsv'
        distance_stats_df.to_csv(total_file, sep='\t', index=False,
                                 columns=['cluster_id', 'chromosome', 'n_peaks',
                                         'mean_dist', 'std_dist', 'compactness', 'cv_dist', 'ratio',
                                         'mean_length', 'cv_length'])
        print(f"  Saved: {total_file} (per-cluster statistics for {len(distance_stats_df)} clusters)")
        saved_files.append(total_file)

    # Save overall statistics
    if save_summary:
        if len(all_distances) > 0:
            overall_mean = np.mean(all_distances)
            overall_std = np.std(all_distances, ddof=1) if len(all_distances) > 1 else 0.0
            overall_cv = overall_std / overall_mean if overall_mean > 0 else 0.0

            summary_df = pd.DataFrame({
                'mean_dist': [overall_mean],
                'std_dist': [overall_std],
                'compactness': [overall_mean],
                'cv': [overall_cv]
            })

            summary_file = output_file.rsplit('.', 1)[0] + '_summary.tsv'
            summary_df.to_csv(summary_file, sep='\t', index=False)
            print(f"  Saved: {summary_file} (overall statistics from {len(all_distances)} pairwise distances)")
            saved_files.append(summary_file)
        else:
            print("  No distances to compute summary statistics")

    return saved_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="3D peak clustering with MM (distance-weighted) centroid calculation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s -i input.narrowPeak -o output.narrowPeak

  # With statistics
  %(prog)s -i input.narrowPeak -o output.narrowPeak --save_total --save_summary

  # Custom max_gap
  %(prog)s -i input.narrowPeak -o output.narrowPeak --max_gap 100

  # Save original peaks by sample with cluster tracking
  %(prog)s -i input.narrowPeak -o output.narrowPeak --save_clustered_peaks

  # Full features
  %(prog)s -i input.narrowPeak -o output.narrowPeak --save_total --save_summary --save_clustered_peaks

Key Features:
  - 3D clustering using (start, summit, end) euclidean distance with max_gap threshold
  - Peak shape (length, asymmetry) naturally incorporated into clustering
  - MM (Majorize-Minimization) distance-weighted centroid calculation
    * Summit: Iterative algorithm where weight = 1/distance (closer points have higher influence)
    * Left/Right boundaries: Same MM distance-weighted approach
  - Optional statistics output with --save_total (per-cluster) and --save_summary (overall)
  - Optional sample-separated peaks with --save_clustered_peaks (includes cluster_id for tracking)
"""
    )
    parser.add_argument("-i", "--input_bed", required=True, help="Input narrowPeak file")
    parser.add_argument("-o", "--output_bed", required=True, help="Output narrowPeak file")
    parser.add_argument("--max_gap", type=int, default=70,
                        help="Maximum distance allowed within a cluster (default: 70)")
    parser.add_argument("--min_cluster_size", type=int, default=1,
                        help="Minimum number of peaks required per cluster (default: 1)")
    parser.add_argument("-@", "--n-threads", type=int, default=16, help="Number of threads (default: 16)")
    parser.add_argument("--save_total", action="store_true",
                        help="Save per-cluster statistics to total.tsv file")
    parser.add_argument("--save_summary", action="store_true",
                        help="Save overall statistics to summary.tsv file")
    parser.add_argument("--save_clustered_peaks", action="store_true",
                        help="Save original peaks by sample to separate files (e.g., output_SAMPLE.narrowPeak)")

    args = parser.parse_args()

    print("="*70)
    print("3D Peak Clustering with MM Distance-Weighted Centroid")
    print("="*70)
    print(f"Input:                  {args.input_bed}")
    print(f"Output:                 {args.output_bed}")
    print(f"Max gap:                {args.max_gap}")
    print(f"Min cluster size:       {args.min_cluster_size}")
    print(f"Threads:                {args.n_threads}")
    print(f"Clustering:             3D euclidean (start, summit, end) + max_gap")
    print(f"Centroid:               MM distance-weighted (weight = 1/distance)")
    print(f"Save total stats:       {args.save_total}")
    print(f"Save summary stats:     {args.save_summary}")
    print(f"Save clustered peaks:   {args.save_clustered_peaks}")
    print("="*70)
    print()

    # Load and validate input
    bed = pd.read_csv(args.input_bed, sep="\t", header=None, dtype={0: str})
    if bed.shape[1] < 10:
        raise ValueError(f"Input must be narrowPeak format (10 columns required, got {bed.shape[1]})")

    # Preprocess data
    bed.columns = [f"col{i}" for i in range(bed.shape[1])]
    bed.rename(columns={
        "col0": "chr", "col1": "start", "col2": "end", "col3": "sample",
        "col4": "raw_score", "col5": "strand", "col6": "signalValue",
        "col7": "pValue", "col8": "qValue", "col9": "peak_offset"
    }, inplace=True)

    bed["pos"] = (bed["start"] + bed["peak_offset"]).astype(int)
    bed["length"] = (bed["end"] - bed["start"]).astype(int)

    bed["raw_score"] = bed["raw_score"].astype(float)
    bed["signalValue"] = bed["signalValue"].astype(float)
    bed["pValue"] = bed["pValue"].astype(float)
    bed["qValue"] = bed["qValue"].astype(float)
    bed["sample"] = bed["sample"].astype(str)

    print(f"Dataset: {bed['chr'].nunique()} chromosomes, {len(bed)} peaks")
    print(f"Processing with {args.n_threads} threads using 3D clustering + MM centroid...")
    print()

    # Process chromosomes in parallel
    chrom_data_list = [
        (chrom, df_chr.reset_index(drop=True), args.max_gap, args.min_cluster_size,
         args.save_total, args.save_summary, args.save_clustered_peaks)
        for chrom, df_chr in bed.groupby("chr")
    ]

    if args.n_threads > 1 and len(chrom_data_list) > 1:
        with Pool(processes=args.n_threads) as pool:
            results = list(tqdm(pool.imap(process_chromosome, chrom_data_list),
                               total=len(chrom_data_list), desc="Chromosomes"))
    else:
        results = [process_chromosome(d) for d in tqdm(chrom_data_list, desc="Chromosomes")]

    # Collect results
    centroid_rows = []
    total_clusters = 0
    valid_clusters = 0
    all_distance_stats = []
    all_pairwise_distances = []
    all_sample_peaks = []

    for centroid_list, n_total, n_valid, stats, distances, sample_peaks in results:
        centroid_rows.extend(centroid_list)
        total_clusters += n_total
        valid_clusters += n_valid
        if args.save_total:
            all_distance_stats.extend(stats)
        if args.save_summary:
            all_pairwise_distances.extend(distances)
        if args.save_clustered_peaks:
            all_sample_peaks.extend(sample_peaks)

    # Save centroid output
    centroid_df = pd.DataFrame(
        centroid_rows,
        columns=["chr", "start", "end", "name", "score", "cluster_id",
                 "signalValue", "pValue", "qValue", "peak"]
    )
    centroid_df = centroid_df.drop_duplicates()
    centroid_df = centroid_df.sort_values(["chr", "start"])
    centroid_df.to_csv(args.output_bed, sep="\t", header=False, index=False)

    # Print results
    print()
    print("="*70)
    print("Results:")
    print(f"  Output:              {args.output_bed}")
    print(f"  Total clusters:      {total_clusters}")
    print(f"  Valid clusters:      {valid_clusters} (>= {args.min_cluster_size} peaks)")
    print(f"  Total centroids:     {len(centroid_df)}")
    print("="*70)

    # Save statistics files
    if args.save_total or args.save_summary:
        print()
        print("Generating statistics files...")
        distance_stats_df = pd.DataFrame(all_distance_stats) if args.save_total and all_distance_stats else None
        saved_files = save_distance_statistics(distance_stats_df, args.output_bed,
                                              all_pairwise_distances,
                                              args.save_total, args.save_summary)

        if saved_files:
            print()
            print("="*70)
            print("Statistics files:")
            for filepath in saved_files:
                print(f"  {filepath}")
            print("="*70)

    # Save sample-specific files
    if args.save_clustered_peaks and len(all_sample_peaks) > 0:
        print()
        print("Saving sample-specific peak files...")
        saved_files = save_clustered_peaks_by_sample(all_sample_peaks, args.output_bed)

        print()
        print("="*70)
        print(f"Sample-specific files saved for {len(saved_files)} samples:")
        for sample, filepath, n_peaks in saved_files:
            print(f"  {sample}: {filepath} ({n_peaks} peaks)")
        print("="*70)

    print("\nDone!")
