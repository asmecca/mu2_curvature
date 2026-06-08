#!/bin/bash

#Nt=$1
declare -a list_Nt=("40" "36" "32" "28" "24" "20" "16")
declare -a list_T=("0.141" "0.156" "0.176" "0.201" "0.235" "0.281" "0.352")
declare -a list_mu=("0.028" "0.056" "0.084" "0.113" "0.14" "0.169" "0.197" "0.225" "0.281")


for j in `seq 0 8`
do    
    for i in `seq 0 6`
    do	
	Nt=${list_Nt[$i]}
	path_to_correlator_directory=/DATA/gen2/${Nt}x24 	
	T=${list_T[$i]}
	mu=${list_mu[$j]}
	echo $Nt
	echo $mu
	echo $T
	if [ -d ${path_to_correlator_directory}/analysis_mu_${mu}/boot ] ; then
	   rm -r ${path_to_correlator_directory}/analysis_mu_${mu}/boot
	fi
	python3 final-analyse-multisource.py ${path_to_correlator_directory} ${Nt} ${mu} ${T}
    done
done
