make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=1
make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=2
DIM is default set to 3 in Make.defs.local
Make.defs.local in this folder is for LLNL LC quartz.
It needs to be placed in Chombo/lib/mk/ 
