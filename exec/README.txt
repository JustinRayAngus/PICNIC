Make options:
make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=1
make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=2
DIM is default set to 3 in Make.defs.local
Make.defs.local in this folder have options for LLNL LC quartz and systems using intel MKL.
It needs to be placed in path_to_chombo/lib/mk/
