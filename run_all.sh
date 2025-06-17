#!/usr/bin/env bash
set -euo pipefail

ratios=(0   0.25 0.5 0.75 1   1.25 1.5 1.75 2   2.25
        2.5 2.75 3   3.25 3.5 3.75 4   4.25 4.5 4.75
        5   5.25 5.5 5.75 6   6.25 6.5 6.75 7   7.25
        7.5 7.75 8   10   15  20   50  75   100)
exe=./main
dataset="Hypergraphs/com-amazon-cmty-hygra"
logfile="Results/result-com-amazon-cmty-hygra.log"

: > "$logfile"

for pct in "${ratios[@]}"; do
  echo "=== $exe $pct $dataset ===" | tee -a "$logfile"
  "$exe" "$pct" "$dataset" 2>&1 | tee -a "$logfile"


  enc_sum=0
  dec_sum=0
  tot_sum=0

  for run in {1..1}; do
    output=$("$exe" "$pct" "$dataset" 2>/dev/null)

    enc_time=$(echo "$output" | grep "Encoding time" | awk '{print $3}')
    dec_time=$(echo "$output" | grep "Decoding time" | awk '{print $3}')
    tot_time=$(echo "$output" | grep "Total running time" | awk '{print $4}')
    

    enc_sum=$(awk "BEGIN {print $enc_sum + $enc_time}")
    dec_sum=$(awk "BEGIN {print $dec_sum + $dec_time}")
    tot_sum=$(awk "BEGIN {print $tot_sum + $tot_time}")
  done

  enc_avg=$(awk "BEGIN {printf \"%.5f\", $enc_sum / 1}")
  dec_avg=$(awk "BEGIN {printf \"%.5f\", $dec_sum / 1}")
  tot_avg=$(awk "BEGIN {printf \"%.5f\", $tot_sum / 1}")

  echo "Average Encoding time: $enc_avg seconds" | tee -a "$logfile"
  echo "Average Decoding time: $dec_avg seconds" | tee -a "$logfile"
  echo "Average Total time   : $tot_avg seconds" | tee -a "$logfile"
  echo "" | tee -a "$logfile"
done
