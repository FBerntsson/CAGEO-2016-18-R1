%Content:
%
% PaperGraphs: This script produces some of the graphs present in the
%    paper. Its intended to illustrate how the functions are intended to 
%    be used.
%
% ThermalModel: This script defines the thermal model used for the 
%    simulations. That is given a set of grid points (X,Y) and a given
%    temperature distribution on the grid it computes the corresponding 
%    thermal conductivity on the grid and also the heat production. This 
%    script should be changed to edit the ThermalModel. It is not used by
%    any of the other functions. All other functions have matrices with the
%    thermal conductivity and heat production as their input.
%
%  Usage:
%     >> [K,A]=ThermalModel(X,Z,T,'on/off')
%
% DirectThermalSolve: This function solves the well-posed direct problem.
%   That is given the heat-flux Qm at the base of the model and the
%   temperature T0 at the surface it computes the temperature distribution
%   T(x,z) for the entire domain. It can also calculate the surface
%   heat-flux. It uses the finite difference method. 
%
%  Usage:
%     >> [T,Q0]=DirectThermalSolve( x , z , HeatCond , HeatProd , T0 , Qm )
%
%
%
% FDM_CreateLinearSystem: The finite difference discretization leads to a
%    linear system of equations A*u=b. This function creates the matrix A. 
%    If changing to a different method of solving the equations
%    this function needs to be edited. It is used by DirectThermalSolve() 
%    and also by LinearTikhonovSolve. 
%
%  Usage:
%     >> [b]=FDM_CreateLinearSystem( x , z , HeatCond )
%
%
% FDM_CreateRighthandside: The finite difference discretization leads to a
%    linear system of equations A*u=b. This function creates the right 
%    handside b. If changing to a different method of solving the equations
%    this function needs to be edited. It is used by DirectThermalSolve().
%
%  Usage:
%     >> [b]=FDM_CreateRighthandside( x , z , HeatCond , HeatProd , T0 , Qm  )
%
%
% LinearTikhonovSolve: Tikhonov regularisation is applied to the inverse 
% geothermal problem and the normal equations are solved using the
% conjugate gradient method. The function only solves the linear problem. 
%
% Usage: 
%  >> [Qtik]=LinearTiknonovSolve( x,z,HeatCond,HeatProd,T0,Q0,Lambda,QtikInit )
%