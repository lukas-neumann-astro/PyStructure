Guide
=======

The configure File
------------------

Here we go through the configuration file step-by-step and explain the necessary user input.

Step 1
^^^^^^
In the first step, the user need to specify some of the key paths as well as the overlay 3D cube to be used.

.. code-block::

    ####################################
    # Step 1: Define the correct Paths #
    ####################################

    # <path to directory with the data files>
    data_dir = "data/"

    # <filename of geometry file>
    geom_file = "List_Files/geometry.txt"
    # <filename of overlay or mask> #should be stored in data_dir
    overlay_file = "_12co21.fits"

    # <Output Directory for Dictionaries>
    out_dic = "Output/"



* ``geom_file`` The geometry file contains information on the source properties (see :ref:`geomFile` for more details).
* ``overlay_file`` The FITS file extension of the overlay 3D data cube to be used. The overlay is used to define the spatial extend of the hexagonal grid.

.. NOTE::
    
    All data filenames must have the following convention:
    
    ``<source_name><file_extension>``
    
    For example: ``ngc5194_co21.fits``, where ngc5194 is the source name and _co21.fits the file extension.
    
    The ``geometry.txt`` file must contain the source name.
    
    
.. _step2:

Step 2
^^^^^^
In the next step, the targeted angular resolution needs to be provided. By default, **this needs to be given in arcsec**.

.. code-block::

    #####################################
    # Step 2: Set the Target Resolution #
    #####################################
    # Set the target resolution for all data in arcseconds (if resolution set to angular)
    target_res = 27.

Step 3
^^^^^^

In this step, you need to define the source or list of sources. Note that all datafiles must start with the source name as specified in this step.

.. code-block::

    #####################################
    # Step 3: Set Sources              #
    #####################################
    #can be single source (e.g. cloud or galaxy), or list of strings: sources = "ngc5194", "ngc5457"
    sources = "ngc5194"

Step 4
^^^^^^

In this step, we define the list of 2D data we want to include in the PyStructure.
It is not required to include any 2D data.

.. code-block::

    #####################################
    # Step 4: Define Bands              #
    #####################################
    # Column 1: short name of band
    # Column 2: description for database
    # Column 3: units
    # Column 4: extension
    # Column 5: path to files
    # Column 6: extension to uc file
    spire250,    SPIRE250,    MJy/sr,    _spire250_gauss21.fits,    ./data/,    _spire250_gauss21_unc.fits
    
.. NOTE::
    
    If no uncertainty file extist, simply write ``not_included`` for the uncertainty file extension.
    
Step 5
^^^^^^

The final step consists of specifying the list of 3D data cubes. Here at least one data cube must be provided.

.. code-block::

    #####################################
    # Step 5: Define Cubes              #
    #####################################
    # Column 1: short name of cube
    # Column 2: description for database
    # Column 3: units
    # Column 4: extension
    # Column 5: path to files
    # Column 6 optional: extension, if 2D map provided
    # Column 7 optional: extension err, if 2D map provided
    12co21, 12CO2-1, K, _12co21.fits, data/
    12co10, 12CO1-0, K, _12co10.fits, data/

.. IMPORTANT::
    
    By default, the **first cube specified in this list will be used as prior line** to define a spectral and spatial mask over which the 3D data cubes will be postprocessed to determine the moment-0 and peak temperature maps. It is therefore adviced to select the brightest line as the first entry.
    
Additional Files
----------------

.. _geomFile:

The geometry file
^^^^^^^^^^^^^^^^^

The ``geometry.txt`` file is a crucial component of the PyStructure. Here, the source properties need to be defined.

.. code-block::

    # column 1: name
    # column 2: R.A. center [decimal degrees J2000]
    # column 3: Dec. center [decimal degrees J2000]
    # column 4: distance [Mpc]
    # column 5: uncertainty in distance
    # column 6: inclination [deg] (LEDA)
    # column 7: uncertainty in inclination
    # column 8: position angle [deg] (LEDA)
    # column 9: uncertainty in position angle
    # column 10: radius 25th magnitude isophote [arcmin]
    # column 11: uncertainty in radius 25th
    # column 12: ref for position angle
    # Attach comments with #
    # Distances and uncertainties from Brent Groves' NED querying program.
    # Inclinations and PA from Maria's compilation
    # Radius from Leda
    #-------------------------------------------------------------
    # !!!Make sure columns separated by only one tab!!!
    #-------------------------------------------------------------
    M82    148.969687    69.679383    3.5    0.2    77.0    NaN    62    NaN    5.60    0.1
    ngc0628    24.1739458    15.7836619    9.0    0.14    7.0    NaN    20.0    NaN    5.00    0.05
    ngc2903    143.042125    21.500833    8.7    0.9    65.0    NaN    204.0    NaN    6.01    0.05
    
    
* The name of the source needs to match the source name of the data files
* The source properties are tailored to match usefull properties characterizing nearby galaxies. Therefore, parameters, such as inclination or position angle, might not be meaningfull if you work on other types of sources (e.g. Galactic clouds). In this case, it is adviced to enter ``0`` or ``NaN`` instead.
