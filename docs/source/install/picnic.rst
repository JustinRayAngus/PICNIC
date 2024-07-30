.. _install-picnic:

.. raw::html

PICNIC
=====

Once you have performed the instrutions in Chombo, Prerequisites, and PETSc (if using PETSc), then follow the instructions here to build PICNIC.

Obtaining PICNIC
----

Checkout PICNIC from LC's czgitlab:

.. code-block:: bash

   git clone ssh://git@czgitlab.llnl.gov:7999/angus1/picnic.git PICNIC

Or, if you don't have an czgitlab account on LC, checkout PICNIC from personal github page:

.. code-block:: bash

   git clone https://github.com/JustinRayAngus/PICNIC.git PICNIC

Make.defs.local
----

Make sure the Make.defs.local file in PICNIC/exec/ folder is placed in Chombo/lib/mk folder and set approporatiely at described in the Chombo installation instructions. 

Build PICNIC
----

Build PICNIC in 1D and in 2D geometries.

.. code-block:: bash

   cd path/to/picnic/exec/
   make -j all MPI=TRUE DEBUG=FALSE DIM=1
   make -j all MPI=TRUE DEBUG=FALSE DIM=2

The various axisymmetric geometries are used by specifying so in the input file at run time.

