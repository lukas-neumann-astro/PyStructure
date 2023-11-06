Installation
============

Getting the scripts
-------------------

The up-to-date version of the PyStructure can be found on the `GitHub page <https://github.com/jdenbrok/PyStructure>`_.
Most easily, you can just clone the repository:

.. code-block:: console

   $ git clone https://github.com/jdenbrok/PyStructure.git



.. _dependancies:

Dependancies
------------

The scripts are based on Python3. See the `requirements.txt` file in the GitHub folder for the list of python modules and their version.
Currently, the list includes:

* astropy==5.2.2
* matplotlib==3.4.3
* numpy==1.22.4
* pandas==2.0.3
* radio_beam==0.3.4
* reproject==0.9.1
* scipy==1.11.3
* spectral_cube==0.6.0


It is recommended to set up a virtual environment. The Python modules can then be installed easily like this:

.. code-block:: console

    (./venv) $ python3 -m pip install -r requirements.txt 

