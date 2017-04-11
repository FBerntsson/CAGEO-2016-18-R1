# CAGEO-2016-18-R1

Title: An efficient regularization method for a large scale ill-posed geothermal problem
Authors: Fredrik Berntsson, Chen Lin, Tao Xu, Dennis Wokiyi
Journal: Computers and Geosciences, 2017.

Contact and Information: Fredrik Berntsson (Email: fredrik.berntsson@liu.se) for questions regadring the Matlab codes.

Abstract:

The inverse geothermal problem consists of finding the temperature and heat-flux in the interior of the earth using 
surface measurements. This is an ill-posed problem that has to be solved using regularization techniques. In the paper
we present an efficient implementation of Tikhonov regularization for solving the inverse geothermal problem. 

The Matlab codes used for creating the graphs in the paper are available in this repository. Read the file Contents.m for 
more information. To get started run PaperGraphs.m which generates most of the figures available in the paper. The codes are
written in such a way that the thermal model can be changes relatively easily. The thermal parameters are seperated from the 
rest of the solution software and the asdsumption is that a function

  >> [K,A]=ThermalModel(X,Z,T)
  
 that calculates matrices K, containing thermal conductivities, and A, containing heat production, for a meshgrid X,Y and 
 temperature T at the grid, is made available. Thus the codes are easily adapted to solve different thermal models.  

Best wishes,

Fredrik Berntsson
