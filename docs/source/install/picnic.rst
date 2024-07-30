.. _install-picnic:

.. raw::html

PICNIC
=====

Once you have performed the instructions in Chombo, Prerequisites, and PETSc (if using PETSc), then follow the instructions here to obtain and build PICNIC.

Obtaining PICNIC
----

If you have an LC czgitlab account, checkout PICNIC from LC's czgitlab:

.. code-block:: bash

   git clone ssh://git@czgitlab.llnl.gov:7999/angus1/picnic.git PICNIC

Or, if you don't have an czgitlab account on LC, checkout PICNIC from github:

.. code-block:: bash

   git clone https://github.com/JustinRayAngus/PICNIC.git PICNIC

Make.defs.local
----

Make sure the Make.defs.local file in PICNIC/exec/ folder is placed in Chombo/lib/mk folder and set approporatiely at described in the Chombo installation instructions. 

Build PICNIC
----

Build PICNIC in 1D and in 2D geometries.

.. code-block:: bash

   cd path/to/PICNIC/exec/
   make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=1
   make -j all MPI=TRUE DEBUG=FALSE OPT=TRUE DIM=2

The various axisymmetric geometries available are used by specifying them in the input file at run time.

