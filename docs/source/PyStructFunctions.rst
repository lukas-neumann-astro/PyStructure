PyStructure Functions
---------------------
The PyStructure class contains a set of functions that help handle and extract the data stored in the dictionary.
You can run the functions using ``database.<function>()``.

.. function:: database.get_coordinates(center = None)

   A function extracting the rightasencion and declination coordinates. If center coordinate is provided, the ra and dec coordimnates are returned as offset in arcsec.

   :param center: optional referece coordinate (e.g., "13:29:52.7 47:11:43"). If provided, the returned values will represent the offset in arccsec with respect to this coordinate.
   :type center: str
   :return: ``ra``, ``dec`` ; two 1D arrays, one for the rightasencion and on for the declination.
   :rtype: np.array
