#!/bin/bash

declare -a list_mu=("0.0" "0.028" "0.056" "0.084" "0.113" "0.140" "0.169" "0.197" "0.225" "0.281")

echo "Do you want to remove previous run? [Yes/No]"
select yn in "Yes" "No"; do
    case $yn in
	Yes ) rm interpolate/boot_pseudo_T.dat; break ;;
	No ) echo "OK"; break ;;
    esac
done
	    
for j in `seq 0 8` #5
do
    mu=${list_mu[$j]}
    echo ${mu}
    python3 final_interpolate.py ${mu} 
done
