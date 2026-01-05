#Modified from https://github.com/MeriamGuellil/aDNA-BAMPlotter (aDNA-BAMPlotter.py v.2.0.1)
#To run:  python aDNA-BAMPlotter_2kb_intervals.py -b in.bam -o out.pdf

import os
import argparse
import datetime
import numpy as np
import pysam
import pysamstats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from tqdm import tqdm
from io import StringIO

parser = argparse.ArgumentParser(description="""Genome-wide BAM Summary Plot (Modified from aDNA-BAMPlotter)""")
parser.add_argument('-b', dest='bamR', required=True, type=str, help='Indexed BAM file (required)')
parser.add_argument('-d', dest='deam', required=False, type=str, help='mapDamage2 misincorporation.txt file (optional)')
parser.add_argument('-o', dest='out', required=True, type=str, help='Output image file (e.g. .pdf, .png)')
parser.add_argument('-q', dest='mqf', required=False, default="30", type=str, help='Mapping quality threshold (default: 30)')
args = parser.parse_args()

date = datetime.datetime.now().strftime("%d/%m/%Y")

#Load BAM
bam = pysam.AlignmentFile(args.bamR, "rb")
head = {n["SN"]:n["LN"] for n in bam.header['SQ']}
#sorted_head = dict(sorted(head.items()))
sorted_head = head #already sorted by ascending node order in bam file


#LOAD MISINCORPORATION
if args.deam:
    deam_df = pd.read_csv(args.deam, header=1, delimiter='\t', encoding='utf-8', skiprows=2)


#LOAD DEPTH
depth = pd.read_csv(StringIO(pysam.depth(args.bamR,"-aa")),delimiter='\t',encoding='utf-8',names=["id","pos","dp"])
depth30 = pd.read_csv(StringIO(pysam.depth(args.bamR,"-aa","-Q",args.mqf)),delimiter='\t',encoding='utf-8',names=["id","pos","dp"])

#create empty containers for later because we need global alignment for a singular plot rather than a plot per contig
offset = 0 #starting the offset at 0 creates a global alignment rather than just contig by contig 
depth_all = []
depth30_all = []
ED_all = []
ED30_all = []

pbar = tqdm(sorted_head.items(), desc="Processing Contigs", unit="contig")

#use offsets to create global alignment 
offsets = {}
for ID, LEN in pbar:
    offsets[ID] = offset

    #Contig data 
    d = depth[depth["id"] == ID].copy()
    d["global_pos"] = d["pos"] + offset
    depth_all.append(d[["global_pos", "dp"]])

    #High qual contigs
    d30 = depth30[depth30["id"] == ID].copy()
    d30["global_pos"] = d30["pos"] + offset
    depth30_all.append(d30[["global_pos", "dp"]])

    #Edit distances (ED)
    NM_all = [r.get_tag("NM") for r in bam.fetch(ID)]
    NM_30 = [r.get_tag("NM") for r in bam.fetch(ID) if r.mapping_quality >= int(args.mqf)]
    ED_all.extend(NM_all)
    ED30_all.extend(NM_30)

    offset += LEN

#Combine coverage
df_cov = pd.concat(depth_all)
df_cov30 = pd.concat(depth30_all)

#bin into 2kb intervals
df_cov["bin"] = (df_cov["global_pos"] // 2000) * 2000
df_cov30["bin"] = (df_cov30["global_pos"] // 2000) * 2000

#Find the median for each 2kb bin
cov_summary = df_cov.groupby("bin")["dp"].median().reset_index(name="MQ0")
cov_summary30 = df_cov30.groupby("bin")["dp"].median().reset_index(name="MQ30")

#Merge both tables
cov_merged = pd.merge(cov_summary, cov_summary30, on="bin", how="outer").fillna(0).sort_values("bin")

#PLOT
plt.rcParams.update({'font.size': 11})
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(20, 6))
fig.suptitle(f"BAM Summary - {os.path.basename(args.bamR)}", fontsize=16, weight='bold')

#Plot 1: Coverage
axs[0].plot(cov_merged["bin"], cov_merged["MQ0"], color="gray", alpha=0.5, label="MQ ≥ 0")
axs[0].plot(cov_merged["bin"], cov_merged["MQ30"], color="#148F77", label=f"MQ ≥ {args.mqf}")
axs[0].fill_between(cov_merged["bin"], cov_merged["MQ30"], alpha=0.3, color="#148F77")
axs[0].set_xlabel("2kb Bins Across Genome")
axs[0].set_ylabel("Median Coverage")
axs[0].legend()
axs[0].set_title("Genome-wide Coverage")
axs[0].set_xlim([cov_merged["bin"].min(), cov_merged["bin"].max()])
axs[0].set_ylim(bottom=0)

#Plot 2: Edit Distance
bins = range(0, max(ED_all + ED30_all) + 1)
axs[1].hist(ED_all, bins=bins, color="gray", alpha=0.5, label="MQ ≥ 0", density=True)
axs[1].hist(ED30_all, bins=bins, color="#148F77", alpha=0.7, label=f"MQ ≥ {args.mqf}", density=True)
axs[1].set_xlabel("Edit Distance")
axs[1].set_ylabel("Frequency")
axs[1].set_title("Edit Distance Distribution")
axs[1].legend()
axs[1].set_xlim([min(bins), max(bins)])
axs[1].set_ylim(bottom=0)

#fix weird extra tick marks below x-axis
from matplotlib.ticker import MaxNLocator, ScalarFormatter

for ax in axs:
    ax.ticklabel_format(style='plain', axis='x')  # turn off sci notation
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))  # clean integer ticks


#Annotate and Save
fig.text(0.99, 0.99, date, ha='right', va='top', fontsize=9)
fig.tight_layout()
fig.savefig(args.out, bbox_inches='tight')
