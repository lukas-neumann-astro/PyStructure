Quick Start
============

Preparing the default
----------------------

The installations of the script comes with a default example which allows you to test whether the code can run on your machine and if you have installed all necessary Python packages.

The first lines of the ``PyStructure.conf`` configure file contain information about the path to the data directory. Adjust the path such that it points to the ``data`` directory that you find within the PyStructure folder:

.. code-block:: console

   ####################################
   # Step 1: Define the correct Paths #
   ####################################

   # <path to directory with the data files>
   data_dir = "<path to data directory>"

   


.. _run_example:

Running the Example
-------------------

Once you have correctly adjusted the ``data_dir`` path directory, you can simply run the script using the following command line: 

.. code-block:: console

   (.venv) $ python3 create_database.py --config PyStructure

The script will print information in the terminal as it runs. The result will be a ``.npy`` file, which contains a Python dictionary with the processed data. It will store this file in the ``Output`` folder (its path is specified in the configure file).

For information on how to open and work with the processed data, see :ref:`Analysis`.