#!/bin/bash

declare -a list_mu=("0.028" "0.056" "0.084" "0.113" "0.140" "0.169" "0.197" "0.225" "0.281")

for j in `seq 0 3`
do
    mu=${list_mu[$j]}
    echo ${mu}
    python3 final_R_ratio.py ${mu}
done
