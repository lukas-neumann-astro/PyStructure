PyStructure Functions
---------------------
The PyStructure class contains a set of functions that help handle and extract the data stored in the dictionary.
You can run the functions using ``database.<function>()``.

.. function:: database.get_coordinates(center = None)

   A function extracting the rightasencion and declination coordinates. If center coordinate is provided, the ra and dec coordimnates are returned as offset in arcsec.

   :param center: optional reference coordinate (e.g., ``"13:29:52.7 47:11:43"``). If provided, the returned values will represent the offset in arccsec with respect to this coordinate.
   :type center: str
   :return: ``ra``, ``dec`` ; two 1D arrays, one for the rightasencion and on for the declination.
   :rtype: np.array


.. function:: database.quickplot_2Dmap(line, s = 50, cmap = None)

    A function that generates a plot showing the integarted intensities of a user defined line.

    :param line: Name of the line in the PyStructure, e.g. "12CO21".
    :type line: str
    :param s: Marker size. User needs to vary if points do overlap or if there is too much space between the scatter points.
    :type s: int
    :param cmap: Colormap, default is "RdYlBu_r"
    :type cmap: int
    :return: 2D scatter plot illustrating the integrated intensities.

.. function:: database.get_vaxis(get_shuff = False)

    A function extracting the spectral velocities

    :param get_shuff: If ``True``, return the shuffled spectral axis.
    :type center: bool
    :return: ``vaxis`` ; 1D arrays with the spectral axis values.
    :rtype: np.array
