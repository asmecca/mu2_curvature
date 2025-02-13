#!/bin/bash

#Nt=$1
declare -a list_Nt=("16" "20" "24" "28" "32" "36" "40")
declare -a list_T=("0.380" "0.304" "0.253" "0.217" "0.190" "0.169" "0.152")
declare -a list_mu=("0.028" "0.056" "0.084" "0.113" "0.14" "0.169" "0.197" "0.225" "0.281")


for j in `seq 0 5`
do    
    for i in `seq 0 6`
    do	
	Nt=${list_Nt[$i]}
	path_to_correlator_directory=/Users/antoniosmecca/Documents/Physics/pdoc_Swansea/mu2/Code/DATA/Gen2L/${Nt}x32 	
	T=${list_T[$i]}
	mu=${list_mu[$j]}
	echo $Nt
	echo $mu
	echo $T
	python3 final-analyse-multisource.py ${path_to_correlator_directory} ${Nt} ${mu} ${T}
    done
done
