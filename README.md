# mu2_curvature
Analysis scripts for the study of the QCD pseudocritical line curvature using lattice QCD described in 2412.20922

# Requirements:
- Any version of Python 3 available
- The 'gvar' library - https://github.com/gplepage/gvar



# Instructions:

1. unzip gen2.zip and gen2l.zip

2. modify 'iter_mu2.sh' in both gen2 and Gen2L directories, inserting the correct path to the directories NtxNs

3. run 'iter_mu2.sh' for both gen2 and Gen2L

4. modify the file 'final_ratio.py' inserting the correct path to the correlators directory NtxNs for both gen2 and Gen2L

5. run 'launch_ratio.sh' for both gen2 and Gen2L

6. run 'launch_interpolate.sh' for both gen2 and Gen2L

7. run 'final_fit_mu.py' to obtain final plot and kappa values

8. run 'plot_interpolate.py' to obtain the plots showing the interpolation of the data as shown in the paper. It takes the value of mu_q as input

The gen2 directory contains also the script 'comp_free_corr.py' to compare the full lattice correlator with the free correlator obtained analytically for the same timeslices as the lattice correlator. Also in this case, one needs to modify the script with the correct path to the gen2 directory NtxNs.

Attention! - ''launch_interpolate.sh'' appends the values of Tpc to a file, so if it is not the first run you should delete the previous run


#OUTPUT:	

In each output directory (analysis_mu_xxx) there will be three subdirectories:
	   	
	- boot : contains .dat files for each channel of interest (vector, axial-vector)
	       	 each .dat file contains all bootstrap samples for all the terms needed for the mu^2 calculation. 

		 The terms are divided in the following columns:

		 0. G(t) (raw, or order 0, correlator) 
		 1. +4 Re < conn.[0] >
		 2. -2 Re < conn.[1] >
		 3. -2 Re < conn.[2] >
		 4. -8.0 < Im(conn.[3]) Im(disc.[0]) >
		 5. +2.0 < conn.[4] * (disc.[1] - disc.[2] + 2*disc.[3]) - <conn.[4]> <disc.[1] - disc.[2] + 2*disc.[3]> >
		 6. Sum of 1., 2., 3., 4., 5.
		 7. 0. + 0.5 * mu^2 * 6. ( G(t) + 0.5 * mu^2* G''(t) )

	- full : contains .dat file for axial-vector and vector channels plus all gamma-matrices options following Luescher definition in openQCD.
	       	 each .dat files contains onlty term 7. from the divisions above averaged over the bootstrap samples. first column is central value second column is standard error.

	- terms: contains .dat file for axial-vector and vector channels plus all gamma-matrices options following Luescher definition.
	  	 each .dat file contains all terms described above averaged over boostrap samples. 
		 There are 16 columns in total, first column is central value of term 0. and second column is standard error of 0. and so on for the next columns