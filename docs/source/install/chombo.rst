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

Chombo
=====

PICNIC uses the Chombo library for data containers and efficient MPI-handling.
Chombo is developed and maintained by LBL-ANAG. Information about Chombo and
instructions for dowloading can be found at <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>

Obtaining Chombo
----

Chombo is controlled using svn. Checking out Chombo requires registering for an anag-repo account at https://anag-repo.lbl.gov. To check out Chombo: 

.. code-block:: bash

   svn co https://anag-repo.lbl.gov/svn/Chombo/trunk/ Chombo

If this is the first time you are checking out Chombo, you may need to checkout using 

.. code-block:: bash

   svn --non-interactive --trust-server-cert --username USERNAME --password PASSWORD co https://anag-repo.lbl.gov/svn/Chombo/trunk Chombo

Specify the path to Chombo
----

There is no need to configure or build Chombo. The libraries needed by PICNIC will be built during the initial compilation of PICNIC. An environment variable CHOMBO_DIR that points to the Chombo directory is expected by PICNIC at compile time:

.. code-block:: bash

   export CHOMBO_DIR="path/to/Chombo/"

Chombo Make.defs.local file
----

Chombo requires a machine-specific Make.defs.local file to be placed in the Chombo/lib/mk/ folder. There is one in the picnic/exec/ folder that works for LC clusters Ruby and Dane at Lawrence Livermore National Laboratory. Copy this file into the Chombo/lib/mk/ directory prior to compiling. For Mac OS and intel MKL, one just needs to change the XTRALIBFLAGS+= line using the appropriate commented out line in that file.

Cleaning Chombo
----

When running into compile issues, it can be necessary to clean the Chombo libraries. This is done using

.. code-block:: bash

   cd path/to/Chombo/lib/
   make realclean
