import gzip
import json

import numpy as np
import pyfastx

def quantile_match_histograms(short_bins, short_counts, long_bins, long_counts):
    """
    Perform quantile (CDF) matching between two length histograms.

    Parameters
    ----------
    short_bins : 1D array
        Unique fragment lengths observed in the SHORT (target) sample,
        sorted in increasing order.
    short_counts : 1D array
        Number of reads/fragments for each length in `short_bins`.
    long_bins : 1D array
        Unique fragment lengths observed in the LONG sample (to be trimmed),
        sorted in increasing order.
    long_counts : 1D array
        Number of reads/fragments for each length in `long_bins`.

    Returns
    -------
    mapped_lengths : 1D array (same shape as `long_bins`)
        For every length value in `long_bins`, the corresponding target
        length (capped so that it never exceeds the original length,
        since we can only trim, not extend reads).
    """
    short_bins = np.asarray(short_bins, dtype=float)
    long_bins  = np.asarray(long_bins, dtype=float)
    short_counts = np.asarray(short_counts, dtype=float)
    long_counts  = np.asarray(long_counts, dtype=float)

    # --- Empirical PMFs & CDFs -------------------------------------------------
    short_pmf = short_counts / short_counts.sum()
    long_pmf  = long_counts  / long_counts.sum()

    short_cdf = np.cumsum(short_pmf)
    long_cdf  = np.cumsum(long_pmf)

    # --- Representative quantile per long bin -----------------------------
    # Use the midpoint between F(x-) and F(x) as the "quantile" of that bin,
    # which corresponds to the average rank of all reads sharing that length.
    long_cdf_prev = np.concatenate(([0.0], long_cdf[:-1]))
    long_quantile  = (long_cdf_prev + long_cdf) / 2.0

    # --- Map quantiles to target lengths via inverse CDF (interpolation) ---
    mapped_lengths = np.interp(long_quantile, short_cdf, short_bins)

    # --- Trim-only constraint: never request lengths > original -----------
    mapped_lengths = np.minimum(mapped_lengths, long_bins)

    return mapped_lengths


# Read fragment length histogram of the ancient sample
with open(snakemake.input.short, "rt") as jsonfile:
    fastp_short = json.load(jsonfile)
total_short = fastp_short['summary']['after_filtering']['total_reads']
hist_short = np.asarray(fastp_short['insert_size']['histogram'])
hist_short = hist_short[29:(np.max(np.where(hist_short > 0)) + 1)]

# Read number of non-merged sequences of the ancient sample
with open(snakemake.input.not_merged, "rt") as f:
    not_merged = int(f.readline().rstrip())

# Read fragment length histogram of the modern sample
with open(snakemake.input.long, "rt") as jsonfile:
    fastp_long = json.load(jsonfile)
total_long = fastp_long['summary']['after_filtering']['total_reads']
hist_long = np.asarray(fastp_long['insert_size']['histogram'])
hist_long = hist_long[29:(np.max(np.where(hist_long > 0)) + 1)]

# Create the rank/quantile match lookup table
mapped_lengths = quantile_match_histograms(np.asarray(range(30, hist_short.shape[0] + 30)),
                                           hist_short,
                                           np.asarray(range(30, hist_long.shape[0] + 30)),
                                           hist_long)
length_lookup = dict(zip(np.asarray(range(30, hist_long.shape[0] + 30)).astype(int), np.round(mapped_lengths).astype(int)))

# Determine the number of non-merged reads and the total reads
prob_not_merged = not_merged / (total_short + not_merged)

# Prepare FastQ files for writing
outfile = gzip.open(snakemake.output[0], "wt")

# Trim FastQs of merged reads
for name, seq, qual, comment in pyfastx.Fastx(snakemake.input.merged, comment=True):
    target_length = length_lookup[len(seq)]
    outfile.write(f"@{name} {comment}\n")
    outfile.write(f"{seq[:target_length]}\n+\n")
    outfile.write(f"{qual[:target_length]}\n")

rng = np.random.default_rng(seed=int(snakemake.wildcards.rep))
# Trim FastQs of non-merged reads
for r1, r2 in zip(pyfastx.Fastx(snakemake.input.pe1, comment=True),
                  pyfastx.Fastx(snakemake.input.pe2, comment=True)):
    if rng.choice([True, False], p=[prob_not_merged, 1-prob_not_merged]):
        outfile.write(f"@{r1[0]} {r1[3]}\n{r1[1]}\n+\n{r1[2]}\n")
        outfile.write(f"@{r2[0]} {r2[3]}\n{r2[1]}\n+\n{r2[2]}\n")
    else:
        target_length = length_lookup[len(r1[1])]
        outfile.write(f"@{r1[0]} {r1[3]}\n")
        outfile.write(f"{r1[1][:target_length]}\n+\n")
        outfile.write(f"{r1[2][:target_length]}\n")
        target_length = length_lookup[len(r2[1])]
        outfile.write(f"@{r2[0]} {r2[3]}\n")
        outfile.write(f"{r2[1][:target_length]}\n+\n")
        outfile.write(f"{r2[2][:target_length]}\n")

outfile.close()
