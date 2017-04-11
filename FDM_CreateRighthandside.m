% FDM_CreateRightHandside: Create the linear system of equations that correspond
% to a finite difference discretization of the direct problem. This
% function only creates the righthandside vector b.  
%
% Usage:
%
% Create the linear system of equations. 
%
function [b]=FDM_CreateRighthandside( x , z , HeatCond , HeatProd , T0 , Qm  )

%
% Find grid parameters
%
 N=length(x);M=length(z);
 dx=x(2)-x(1);dz=z(2)-z(1);
%
% Create righthandside.
%
 b=zeros(N*M,1);
%
% First interior points. The point we discretize is (x(i),z(j)) and, e.g.
% k(x(i),y(j))=K(j,i), for the matrix created by Thermal model
%
 for i=2:N-1,
   for j=2:M-1,       
      k=i+(j-1)*N;
      %
      % Compute coefficient values at half-points using linear
      % interpolation
      %
       b(k)=-HeatProd(j,i);
   end;
 end;


%
% At the surface z=0 we use temperature boundary conditions
% 
 j=1;
 for i=1:N,
     k=i+(j-1)*N;
     b(k)=T0(i);
 end;
%
% At the base of the model we have a specified heat flux
 j=M;
 for i=1:N,
     k=i+(j-1)*N;
     b(k)=Qm(i)/HeatCond(j,i); 
 end;  

end