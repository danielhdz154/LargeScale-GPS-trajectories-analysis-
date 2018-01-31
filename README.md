# LargeScale-GPS-trajectories-analysis-

This repository contains all the code used for processing large-scale GPS trajectory data (crowdsourcing navigation app data) as part of the research subject Infrastructure Engineering Research Project CVEN900064 for the Master of Engineering (Spatial) at Melbourne University. The output of the code was used for the report 'Evaluating the usability of large-scale datasets as ground truth for urban change and traffic incident detection'. Contact information: hernandezd@student.unimelb.edu.au

The following notebooks demonstrate the workflow used for processing and evaluating the trajectory data for a study case in Melbourne, Australia using raw trajectory data from 2016: 

> 001 - Read and filter GPS data.ipynb <br>
> 002 - Map matching process.ipynb <br>
> 003 - Join map matched data .ipynb <br>
> 004 - Determine traffic normal conditions for peak times.ipynb <br> 
> 005 - Time series modelling - Detecting potential changes in road network.ipynb <br>

Required python scripts, packages and libraries:

-  mapmatcher.py - Provided in the repository
-  pandas
-  numpy 
-  sqlalchemy
-  GDAL 
-  fiona
-  shapely
-  networkx
-  rtree
-  matplotlib

