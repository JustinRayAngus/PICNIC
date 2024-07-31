.. PICNIC documentation master file, created by
   sphinx-quickstart on Sun Jul 28 10:00:30 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:orphan:

PICNIC
-----

PICNIC is a framework for implementing and using advanced **Particle-In-Cell** and **Monte-Carlo collision** methods.

The following electromagnetic PIC time solvers are implemented:

    - explicit leap frog
    - energy-conserving semi-implicit
    - energy-conserving fully implicit

It has modules for the following collisional processes:

    - Coulomb collision
    - charge-exchange
    - elastic collisions between neutrals (hard sphere and variable hard sphere)
    - null method for electron tranport in a gas (elastic + excitation + ionization)
    - nuclear fusion

The code supports the following geometries:

    - 1D/2D planar
    - 1D/2D axisymmetric cylindrical
    - 1D axisymmetric spherical

.. toctree::
   :caption: INSTALLATION
   :maxdepth: 1
   :hidden:

   install/prerequisites
   install/chombo
   install/petsc
   install/picnic

.. toctree::
   :caption: Run PICNIC
   :maxdepth: 1
   :hidden:

   runpicnic/onlc

