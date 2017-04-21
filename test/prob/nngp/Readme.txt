# Introduction:

This is a folder for testing the log-likelihood of Nearest 
Neighbor Gaussian Process(NNGP).


# How to run the test:

The makefile in this folder is written in a relative path. One can
execute the test program in a command window under this directory. 
Use “make” to compile the object file “nngp_test.out” and then use 
command “./nngp_test.out” to get the result.


# Folder Structure:

data: 		HMCprep.R: 		Generating the Synthetic data and other input
					including index and distance matrix of the 						nearest neighbors
		check-nnIndx.R:		Checking the generated nearest index matrix 

lnngp_test: 	nngp_log_test.cpp: 	Test file of nngp-log.hpp
		check_lnngp.R		Compare the results with code R


 






