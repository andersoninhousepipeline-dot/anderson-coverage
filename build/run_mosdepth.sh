#!/usr/bin/env bash
# Run mosdepth over the reference panel for every reference BAM.
set -euo pipefail
MOS=/home/nextflowserver/miniconda3/envs/mosdepth_env/bin/mosdepth
BED=/data/bed/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated.bed
BAMDIR=/data/bed/reference
OUT=/data/bed/webapp/build/mosdepth
THR=1,10,20,30,50,100
MAXJOBS=9
cd "$OUT"
run_one() {
  local slug="$1" bam="$2" base
  base=$(basename "$bam" .bam)
  if [ -f "${slug}__${base}.thresholds.bed.gz" ]; then echo "skip $slug/$base"; return; fi
  "$MOS" -t 4 --by "$BED" --thresholds "$THR" --no-per-base --fast-mode \
    "${slug}__${base}" "$bam" 2>/dev/null
  echo "done $slug/$base"
}
declare -A DIRSLUG=( [Normal-Male]=normal-male [Normal-Female]=normal-female
  [Male-Inf]=male-infertility [Female-Inf]=female-infertility [AF]=af [POC]=poc )
for d in "${!DIRSLUG[@]}"; do
  for bam in "$BAMDIR/$d"/*.bam; do
    while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done
    run_one "${DIRSLUG[$d]}" "$bam" &
  done
done
wait
echo "ALL_MOSDEPTH_DONE"
