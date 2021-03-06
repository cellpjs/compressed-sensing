Implements a compressed sensing algorithm
=========================================

What is compressed sensing?
---------------------------

Consider the matrix equation  
r = Ax  
where  
x is a M x 1 signal vector which is unknown but assumed to be sparse (small number of nonzero elements)  
A is a N x M measurement matrix with N < M (we get to design this matrix)  
r is a N x 1 measurement vector (it is a projection of x into a smaller dimension, thus compressed)  
The goal is to reconstruct x from r

What kind of measurement matrix A do we use?
--------------------------------------------

We use a Low Density measurement matrix with a few fixed number of 1's in each column and each row, and with the remaining elements equal to 0. This enables us to use an efficient reconstruction algorithm.

How to test the algorithm
-------------------------

main_test.m is the main matlab script for testing the algorithm. It does the following:  
1. Call the function generate_ldf() defined in ldf_generate/generate_ldf.m which generates a low density measurement matrix  
2. Randomly generate sparse signal x and measurement r = Ax + n where n is measurement noise  
3. Calculate performance bound where a genie tells you where the nonzero elements of x are  
4. Run the reconstruction algorithm by calling suprem_decoder() defined in suprem_decoder.m which is a wrapper function for the actual decoder written in C (suprem.c)  
5. Compare the performance of the algorithm against the genie-aided one  

Note: Before running suprem_decoder, while MATLAB is in the subfolder called decoder, enter the following command:  
`mex suprem.c`

REFERENCE
---------
A Coding Theory Approach to Noisy Compressive Sensing Using Low Density Frames  
Mehmet Akçakaya, Jinsoo Park, and Vahid Tarokh  
IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 59, NO. 11, NOVEMBER 2011  
freely available version at  
http://arxiv.org/pdf/0903.0650.pdf