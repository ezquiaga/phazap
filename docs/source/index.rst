.. phazap documentation master file, created by
   sphinx-quickstart on Wed Nov  1 06:35:36 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to phazap's documentation!
==================================

**Phazap** is a package to efficiently post-process gravitational-wave (GW)
parameter estimation data to obtain the phases and polarization state of
the signal at a given detector and frequency. Details on the method
are presented in `Ezquiaga, Hu, Lo (2023) <https://arxiv.org/abs/2308.06616>`_. 

Installation
------------
The package can be installed using pip:

.. code-block:: bash

  pip install phazap


or directly from the github repository for the latest features and updates:

.. code-block:: bash

  git clone https://github.com/ezquiaga/phazap.git
  cd phazap
  pip install .

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   examples
   main_api
   full_api

