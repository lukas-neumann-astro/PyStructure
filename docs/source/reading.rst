.. _Analysis:

Working with PyStructures
=========================

How to open the file
--------------------

When you run the ``create_database.py`` script (see :ref:`run_example`), it creates a
``.npy`` file, which simply contains a dictionary. To open this dicionary (for example in another ``Python`` scripts
or in Jupyter Notebook), you can simply run:

.. code-block::

  import numpy as np
  database = np.load("<path to npy file>", allow_pickle = True).item()

The instance ``database`` is now a ``Python`` dictionary, and the relevant, processed data can be extracted using this infrastructure.
For example ``database['rgal_kpc']`` returns an array containing the galacocentric distances of each point (in kpc).

The PyStructure class
---------------------
The download comes with a script that lets you handle the PyStructure output as a ``Python`` class.
This makes it easier to work and use with the databases.

To read in a database, you can use the following synthax:
.. code-block::

  import sys
  sys.path.append("<path to PyStructure/scritps folder")
  import PyStructure as ps

  database = ps.PyStructure("path_to_file.npy")

This way you have extracted the databse as a dictionary, and can access it (for example the galactocentric radii) using ``database.struct['rgal_kpc']``.

Quick Examples
--------------

Using the PyStructure class environment, some useful commands are:

* **List of lines (3D cubes) included in PyStructure**
.. code-block::

  print(database.lines)
  >>> ['12CO21', '12CO10']

* **Extract and plot spectrum (e.g. of brightest sightline)**
.. code-block::

  import matplotlib.pyplot as plt

  vaxis = database.get_vaxis()

  #extract the CO(1-0) integrated intensities
  ii_co10=database.struct['INT_VAL_12CO10']

  #find index of largest intensity (i.e. the brightest sightline)
  idx_brightes = np.argmax(ii_co10)

  spec_co10 = database.struct['SPEC_VAL_12CO10'][idx_brightes,:]

  #plot spectrum
  plt.figure()
  plt.step(vaxis, spec_co10)
  plt.show()

.. image:: spec.png
     :width: 600

* **Make 2D map of integrated intensities (quicklook)**

.. code-block::

  database.quickplot_2Dmap('12CO10')

* **Make 2D map of integrated intensities (more extended)**

  .. code-block::

    database.quickplot_2Dmap('12CO10')

PyStructure Functions
---------------------

REGISTER_VARIABLE("example_variable_documentation(param1, param2, ...)",
                  functionNameHere, R"DOC(
Example variable documentation showing the use of markup.
All the reStructuredText commands as in python function docstrings are valid
here and since it is a raw string we can wrap the string however we see fit:

* we can also add bullet lists
  with multi line entries and *inline* **markup**
* we can add links to other variables and formulas: :math:`x^2`
* and provide documentation for the parameters

Warning:
  Also warnings or info messages

See Also:
  :b2:var:`example_variable_documentation` (linking to ourselves :D)

Parameters:
  param1 (int): parameter description
  param2 (str): and even more parameter descriptions
    across various lines
    and even in **bold**

    Empty lines start a new paragraph in here

)DOC");
