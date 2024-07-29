.. _install-chombo:

.. raw::html

   <style>
   .rst-content section>img {
       width: 30px;
       margin-bottom: 0;
       margin-top: 0;
       margin-right: 15px;
       margin-left: 15px;
       float: left;
   }
   </style>

.. only:: html

   .. image:: chombo_image.jpeg

Obtaining Chombo
=====

PICNIC uses the Chombo library for data containers and efficient MPI-handling.
Chombo is developed and maintained by LBL-ANAG. Information about Chombo and
instructions for dowloading can be found at <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>

Chombo is controlled using svn. Checking out chombo require registering for an anag-repo account at https://anag-repo.lbl.gov. To check out Chombo: 

.. code-block:: bash

   svn co https://anag-repo.lbl.gov/svn/Chombo/trunk/ chombo

If this is the first time you are checking out chombo, you may need to checkout using 

.. code-block:: bash

   svn --non-interactive --trust-server-cert --username USERNAME --password PASSWORD co https://anag-repo.lbl.gov/svn/Chombo/trunk chombo

Specify the path to chombo
----

There is no need to configure or build Chombo. The libraries needed by PICNIC will be built during the initial compilation of PICNIC. An environment variable CHOMBO_DIR that points to the chombo directory is exected by PICNIC at compile time. As an alternative, one may manually specify the path to chombo/lib in the exec/GNUmakefile.

Chombo Make.defs.local file
----

Chombo requires a machine-specific Make.defs.local file to be placed in the chombo/lib/mk/ folder. There is one in the picnic/exec/ folder that works for LC clusters Ruby and Dane at Lawrence Livermore National Laboratory. Copy this file into chombo/lib/mk/ prior to compiling.

Cleaning Chombo
----

When running into compile issues, it is often necessary to clean the Chombo libraries. This is done using

.. code-block:: bash

   cd path/to/chombo/lib/
   make distclean
