.. _Analysis:

Working with PyStructures
=========================

How to open the file
--------------------

When you run the ``create_database.py`` script (see :ref:`run_example`), it creates a
``.npy`` file, which simply contains a dictionary. To open this dicionary (for example in another ``Python~~ scripts
or in Jupyter Notebook), you can simply run:

.. code-block::

  import numpy as np
  database = np.load("<path to npy file>"", allow_pickle = True).item()

The instance ``database`` is now a ``Python`` dictionary, and the relevant, processed data can be extracted using this infrastructure.
For example ``database['rgal_kpc']`` returns an array containing the galacocentric distances of each point (in kpc).

The PyStructure class
---------------------
