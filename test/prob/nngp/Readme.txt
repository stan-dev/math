# Introduction:

This is a folder for developing code of the log-likelihood of Nearest 
Neighbor Gaussian Process(NNGP).


# How to run the test:

The makefile in this folder is written in a relative path. To execute the 
test program, change the directory into the folder “nngp_dev” in the 
terminal and use “make” to compile the object file “a.out”. 
Use command “./a.out” to get the result.


# Folder Structure:

data: 		HMCprep.R: 		Generating the Synthetic data and other input
					including index and distance matrix of the 						nearest neighbors
		check-nnIndx.R:		Checking the generated nearest index matrix 

lnngp:		nngp-log.hpp:		Header file of calculating the log-likelihood 					of NNGP

lnngp_test: 	nngp_log_test.cpp: 	Test file of nngp-log.hpp
		check_lnngp.R		Compare the results with code R


 






