%LinearTikhonovSolve: Tikhonov regularisation is applied to the inverse 
% geothermal problem and the normal equations are solved using the
% conjugate gradient method. 
%
% Usage: 
%  >> [Qtik]=LinearTiknonovSolve( x,z,HeatCond,HeatProd,T0,Q0,Lambda,QtikInit )
%
%
function [Qtik]=LinearTikhonovSolve( x,z,HeatCond,HeatProd,T0,Q0,Lambda,QtikInit )

%
% Split the solution in two parts T=T1+T2, where T1 is computed first, 
% and the equation for T2 is linear. 
%
 [T1,Q1]=DirectThermalSolve( x , z , HeatCond , HeatProd , T0 , zeros(size(Q0)) );

%
% Now we solve the linear ill-posed equation K*Qm=Q0-Q1 using Tihkonov
% regularization. First create a create the linear system of equations 
% and compute LU-decomposition. Computing the LU decomposition here is 
% a relatively large investment but later each CG step becomes cheaper
% as the decomposition is reused.
%
 [A]=FDM_CreateLinearSystem( x , z , HeatCond ); 
 [L,U,P]=lu(A);

%
% Create a matrix such that b=W*Qm is the right-hand-side of the finite
% difference equations.
%
 N=length(x);M=length(z);dz=z(2)-z(1);
 W1=[sparse((M-1)*N,N);sparse(diag( HeatCond(M,:).^-1)) ];

%
% Now create a matrix such that Q0=W2*T
% 
 W2=sparse(diag(HeatCond(1,:)))*[speye(N),-speye(N),sparse(N,(M-2)*N)]/dz;

%
% Now we have K=W2*inv(A)*W1*Qm and K^T=W1^T*inv(A^T)*W2^T. Both can be 
% evaluated efficiently using the LU decomposition. Proceed to evaluate 
% K^T*(Q0-Q1). This is the right-handside in the normal equations.
%
 b=W1'*(P'*(L'\(U'\(W2'*(Q0-Q1)))));
 
%
% Now create an anonymous function handle that corresponds to
% (K^T*K+lambda^2*I). This is used for the CG solution. Note that here
% we use the LU decomposition but if the matrix A is too large we can 
% switch to using gmres() instead.
% 
 tikfun=@(Qk) (W1'*(P'*(L'\(U'\(W2'*(W2*(U\(L\(P*(W1*Qk)))))))))+Lambda^2*Qk);
 
%
% Finally solve the normal equations using the Conjugate gradient method.
% If we have a better initial guess then use that instead of zero. The
% stopping criteria is pretty harse and can be lowered a fair bit.
%
 if nargin<8,QtikInit=zeros(size(Q0));,end
 Qtik = pcg(tikfun,b,1e-8,50,[],[],QtikInit);
 
end
