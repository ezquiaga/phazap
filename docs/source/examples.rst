Examples
========

Example 1: a type-I and a type-II lensed GW signals from the same binary black hole (BBH) merger
------------------------------------------------------------------------------------------------
.. note::
  Before running this example, please make sure that
  you have pulled the data files

  * examples/event_1_PE_samples.json
  * examples/event_2_PE_samples.hdf5

  from our github repository using ``git lfs``. Once you have ``git-lfs`` installed,
  you can do this by running

  .. code-block:: bash

    $ git lfs pull

  in the root directory of the repository.

.. include:: tutorial_example_1.rst

Example 2: analyzing GW191215_223052 and GW191222_033537 from GWTC-3 public data release
----------------------------------------------------------------------------------------
.. note::
  Before running this example, please make sure that
  you have downloaded the PE data release files associated to 
  these two events from `Zenodo <https://zenodo.org/records/8177023>`_.

  Alternatively you can download them from the terminal using `wget`

  .. code-block:: bash

    $ wget https://zenodo.org/records/8177023/files/IGWN-GWTC3p0-v2-GW191215_223052_PEDataRelease_mixed_cosmo.h5
    $ wget https://zenodo.org/records/8177023/files/IGWN-GWTC3p0-v2-GW191222_033537_PEDataRelease_mixed_cosmo.h5

.. include:: tutorial_example_2.rst