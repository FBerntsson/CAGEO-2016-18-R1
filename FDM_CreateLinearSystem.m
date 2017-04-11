% FDM_CreateLinearSystem: Create the linear system of equations that correspond
% to a finite difference discretization of the direct problem. This
% function only creates the matrix and not the righthandside. 
%
% Usage:
%  >>[A]=FDM_CreateLinearSystem( x , z , HeatCond , HeatProd , T0 , Qm  );
%
function [A]=FDM_CreateLinearSystem( x,z,HeatCond ); %, HeatProd , T0 , Qm  )

%
% Find grid parameters
%
 N=length(x);M=length(z);
 dx=x(2)-x(1);dz=z(2)-z(1);

%
% Create matrix. The most efficient way to build the matrix in Matlab
% is to create vectors with row index, column index and the non-zero 
% elements. We have roughly 5 non-zero elements in each row so the 
% size of the vectors is predictable. A pointer ptr is used to indicate
% where the next non-zero element is to be stored in the vectors.
%

% A=sparse(N*M,N*M);

 Elements=zeros(5*M*N,1);Rows=zeros(5*M*N,1);Columns=zeros(5*M*N,1);
 ptr=1;

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
       Coefs = ([HeatCond(j-1,i) HeatCond(j,i-1) HeatCond(j,i+1) HeatCond(j+1,i)]+HeatCond(j,i))/2;       
       Rows(ptr:ptr+4)=[k k k k k]';
       Columns(ptr:ptr+4)=[k-N k-1 k k+1 k+N]';
       Elements(ptr:ptr+4)=[ Coefs(1)/dz^2,Coefs(2)/dx^2,-(Coefs(1)+Coefs(4))/dz^2-(Coefs(2)+Coefs(3))/dx^2,Coefs(3)/dx^2,Coefs(4)/dz^2]';
       ptr=ptr+5;
   end
 end

%
% Fix the derivative boundary conditions at x=0 and x=L1
%
 for j=2:M-1;
     k1=1+(j-1)*N;
     kN=N+(j-1)*N;
     Coefs=[1,-1];
     Rows(ptr:ptr+3)=[k1 k1 kN kN]';
     Columns(ptr:ptr+3)=[k1 k1+1 kN kN-1]';
     Elements(ptr:ptr+3)=[Coefs/dx Coefs/dx]';
     ptr=ptr+4;
 end;

%
% At the surface z=0 we use temperature boundary conditions
% 
 j=1;
 for i=1:N,
     k=i+(j-1)*N;
     Rows(ptr)=[k]';
     Columns(ptr)=[k]';
     Elements(ptr)=[1]';
     ptr=ptr+1;  
 end;

%
% At the base of the model we have a specified heat flux
%
j=M;
 for i=1:N,
     k=i+(j-1)*N;
     Coefs=[1,-1];
     Rows(ptr:ptr+1)=[k k]';
     Columns(ptr:ptr+1)=[k k-N]';
     Elements(ptr:ptr+1)=[Coefs/dz]';
     ptr=ptr+2;
 end;  
 
 %
 % Now create the sparse matrix from the vector of non-zero elements.
 %
 A=sparse( Rows(1:ptr-1),Columns(1:ptr-1),Elements(1:ptr-1));

end
