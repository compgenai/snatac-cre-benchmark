#!/usr/bin/env python3
"""
≥50% reciprocal overlap dedup (CD-HIT-style greedy).

Algorithm:
  1. Sort peaks by MACS2 -log10 pValue (narrowPeak col 8) desc.
     [Note] narrowPeak col 7 is signalValue and col 8 is -log10 pValue;
     the sort key `-k8,8gr` therefore ranks by -log10 pValue, not signalValue.
  2. Assign rank 1..N.
  3. bedtools intersect -a ranked -b ranked -wa -wb -f 0.5 -r
     -> pairs with reciprocal 50% overlap (includes self-pairs).
  4. Filter out self-pairs (rank_a == rank_b).
  5. Build adjacency: rank -> list of overlapping ranks.
  6. Greedy: iterate by rank, keep peak if not already marked redundant,
     mark all higher-rank overlapping peaks as redundant.
  7. Output kept peaks (original narrowPeak format).

Usage:
  python dedup_50pct.py <input.narrowPeak> <output.narrowPeak> <tmpdir>
"""
import os, sys, time, subprocess
from pathlib import Path
from collections import defaultdict

if len(sys.argv) != 4:
    print("Usage: dedup_50pct.py <input.narrowPeak> <output.narrowPeak> <tmpdir>")
    sys.exit(1)

IN  = Path(sys.argv[1])
OUT = Path(sys.argv[2])
TMP = Path(sys.argv[3]); TMP.mkdir(parents=True, exist_ok=True)

BEDTOOLS = os.environ.get("BEDTOOLS", "bedtools")   # override: export BEDTOOLS=/path/to/bedtools

t0 = time.time()
def log(msg): print(f"[{time.time()-t0:6.1f}s] {msg}", flush=True)

# --- 1) Sort by col 8 (signalValue) desc, assign rank ---
ranked_full = TMP / f"{IN.stem}.ranked.narrowPeak"  # 11-col (10 + rank)
ranked_bed  = TMP / f"{IN.stem}.ranked.bed"          # 4-col (chr,start,end,rank) for intersect

log(f"[1] sorting by signalValue desc: {IN}")
subprocess.run(
    f"sort -k8,8gr {IN} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $0, NR}}' > {ranked_full}",
    shell=True, check=True,
)
subprocess.run(
    f"awk 'BEGIN{{OFS=\"\\t\"}}{{print $1, $2, $3, $11}}' {ranked_full} | "
    f"sort -k1,1 -k2,2n > {ranked_bed}",
    shell=True, check=True,
)
N = int(subprocess.check_output(["wc", "-l", str(ranked_full)]).split()[0])
log(f"    total peaks: {N}")

# --- 2) self-intersect with reciprocal 50% overlap, skipping self pairs ---
pairs = TMP / f"{IN.stem}.pairs.tsv"
log("[2] bedtools intersect -f 0.5 -r (self-self)")
# include -sorted for efficiency
subprocess.run(
    f"{BEDTOOLS} intersect -a {ranked_bed} -b {ranked_bed} -wa -wb -f 0.5 -r -sorted "
    f"| awk '$4 != $8 {{print $4, $8}}' > {pairs}",
    shell=True, check=True,
)
log(f"    pair file built: {pairs} ({pairs.stat().st_size/1e6:.1f} MB)")

# --- 3) build adjacency list ---
log("[3] building adjacency")
adj = defaultdict(list)
with open(pairs) as f:
    for line in f:
        a, b = line.split()
        a, b = int(a), int(b)
        # only store directed edges from smaller rank -> larger (to save memory)
        if a < b:
            adj[a].append(b)
        else:
            adj[b].append(a)
log(f"    {sum(len(v) for v in adj.values())} directed edges across {len(adj)} unique peaks")

# --- 4) greedy: iterate by rank, mark redundant ---
log("[4] greedy dedup")
redundant = bytearray(N + 1)  # rank in 1..N, fast bool array
n_kept = 0
for rank in range(1, N + 1):
    if redundant[rank]:
        continue
    n_kept += 1
    for other in adj.get(rank, ()):
        if other > rank:
            redundant[other] = 1
log(f"    kept: {n_kept} / {N}  ({100*n_kept/N:.1f}%)")

# --- 5) write output: only kept peaks (using original narrowPeak format, drop rank col) ---
log(f"[5] writing output: {OUT}")
with open(ranked_full) as fin, open(OUT, "w") as fout:
    for i, line in enumerate(fin, start=1):
        if not redundant[i]:
            # drop last (rank) column
            parts = line.rstrip("\n").split("\t")
            fout.write("\t".join(parts[:-1]) + "\n")
# verify
n_out = int(subprocess.check_output(["wc", "-l", str(OUT)]).split()[0])
assert n_out == n_kept, f"mismatch: kept={n_kept}, written={n_out}"

# --- cleanup intermediate ---
for f in (ranked_full, ranked_bed, pairs):
    try: f.unlink()
    except OSError: pass

log(f"[DONE] {IN.name}: {N} -> {n_kept} ({100*n_kept/N:.1f}% kept, {100*(N-n_kept)/N:.1f}% removed)")
