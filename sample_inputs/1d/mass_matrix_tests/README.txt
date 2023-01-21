These decks are used for verification of the mass matrices.

Need to make realclean and then compile with the
MASS_MATRIX_TEST flag defined in the GNUmakefile

The code will run JFNK using particles for J, but it will also compute the 
mass matrices and compute a test J using them that is written to the output file. 
A direct comparison of these J calculations can be made. See picnic_massMatrixTesting_1D.m
and picnic_massMatrixTesting_2D.m

