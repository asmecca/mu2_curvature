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

8. run 'plot_interpolate.py' to obtain the plots showing the interpolation of the data as shown in the paper

Attention! - ''launch_interpolate.sh'' appends the values of Tpc to a file, so if it is not the first run you should delete the previous run