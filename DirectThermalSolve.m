%DirectThermalSolv: Solve the direct heat conduction problem for the 
% given coefficients and boundary data. The finite difference equations
% are created by two subroutines that has to be provided. 
%
% Usage:
%  [T,Q0]=DirectThermalSolve( x , z , HeatCond , HeatProd , T0 , Qm );
% 
function [T,Q0]=DirectThermalSolve( x , z , HeatCond , HeatProd , T0 , Qm );
 
%
% Set grid
%
 Lx=max(x);N=length(x);dx=x(2)-x(1);
 Lz=max(z);M=length(z);dz=z(2)-z(1);
 
%
% Create the linear system of equations and the right hand side.
%
 [A]=FDM_CreateLinearSystem( x,z,HeatCond );
 [b]=FDM_CreateRighthandside( x , z , HeatCond , HeatProd , T0 , Qm );

%
% Solve the finite difference equations. Here there are several options.
% For small matrices direct methods work well. Generally the LU
% decomposition but if the matrix is designed to be symmetric (which should
% possible for a symmetric finite difference stencil) then use the Cholesky
% decomposition instead. If the grid size is too large then instead use a
% Krylov subspace method, e.g. gmres() or pcg().
%

 T=reshape(A\b,N,M)'; 

%
% Calculate heat-flux at surface level by a first order accurate 
% difference quotient.
%
 Q0=(HeatCond(1,:).*(T(1,:)-T(2,:)))'/dz;
 

end
