.. _install-chombo:

.. raw::html

   <style>
   .rst-content section>img {
       width: 20px;
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
