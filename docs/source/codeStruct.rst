Code Structure
==============

The PyStructure File structure
------------------------------

When you clone or install the GitHub repository, you need to keep the following files and folders: 

| PyStructure
| ├── ListFiles 
| │   ├── band_list.txt
| │   ├── cube_list.txt
| │   └── geometry.txt
| ├── scripts
| │   └── <all scripts in here>
| ├── PyStructure.conf
| └── create_database.py


The key files here that require user input to customize to your project are:

* ``geometry.txt``: This file contains key information about the sources in your project

* ``PyStructure.conf``: The configuration file. Here you need to specify the key details of your dataset to make the script work


How does the PyStructure work?
------------------------------

The script can be broken down into two main steps:

#. Homogenize the Data:
	Based on a user-defined field-of-view and angular resolution, the other datasets will be reprojected and convolved to match sightline-by-sightline. The resulting construct will be a table with a list of sightlines with corresponding coordinates, line spectra or 2D map intensities.

#. Process the 3D Data:
	The code will process the 3D data cubes. This includes:

	* Determining a signal masked based on a user-defined prior line; 
	* Computing S/N-optimized moment 0 (intensity), 1 (mean velocity), 2 (line width), and 8 (peak intensity) maps;
	* Shuffling the spectra by the line-of-sight velocity, such that stacking can be performed.

Intension behind the Code
-------------------------
This code framework has been developed for dealing with a set of spectroscopic 3D astrophysical data spanning different frequencies. The challenging in analyzing such a dataset consists of adequately homogenizing them such that a pixel-by-pixel investigation is possible. Furthermore, the code is optimized to compute key products of 3D data cubes to make a robust analysis accessible. Cubes and maps can be combined.