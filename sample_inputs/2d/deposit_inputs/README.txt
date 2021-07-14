sample inputs for testing the deposit method for depositing
the charge and current of a particle to different locations on
the grid. Only the cloud-in-cell (CIC) method is implemented
correctly for depositing to all grid location (cell center, 
face center, edge center, and node center)

test0_2D - A single particle moving diagonally at a fixed speed
           across a 2D periodic domain. CIC deposit method is used
           and outputs are given for all data container types.
           The results can be used to verify by hand the deposit
           method (CIC) is correct and that MPI is working

test1_2D - Same as test0_2D, but with 4 particles instead of 1 

 
