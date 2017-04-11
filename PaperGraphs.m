%=========================================================================
%
% First setup the Numerically generated test problem. The direct 
% problem is solved iteratively until the temperagture distribution 
% has converged. We print the graphs for the subsection about the 
% test problem. Also after the code has finished we have surface 
% heat-flux Q0 to use as input for the inverse problem solvers.
%
%=========================================================================


% Set up the domain and grid 

 Lx=400e3;         % Length of the model [m] 
 Lz=80e3;          % Deapth of the model [m]

 N=300;x=Lx*(0:N-1)/(N-1); % Grid in x and z.
 M=250;z=Lz*(0:M-1)/(M-1);
 [X,Z]=meshgrid(x,z);

%
% Set boundary conditions for the numerical model.
%

 Qm=(60-15*cos(2*pi*x/Lx)+8*sin(3*pi*x/Lx)-3*sin(5*pi*x/Lx))*1e-3/1.9; 
 T0=10*ones(size(x));  
 plot(x/10^3,10^3*Qm,'LineWidth',1.5);
 xlabel('Horizontal Coordinate: x [ km ]','FontSize',14);
 ylabel('Heat-flux: Q_m [ mW/m ]','FontSize',14)

 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'Boundary conditions for the numerical example.\n\n');
 input('Press return to continue');


%=====================================================================% 
%
% First solve the direct problem by fixed point iteration. This means
% that we assume an initial guess for the temperature distribution 
% T(x,z) and compute the corresponding values for heat production Ap
% and thermal conductivity kappa. In each step we solve a linear problem
% for the new guess for the temperature and recompute the heat production
% and thermal conductivity. 
%
% After convergence we have the solution to the non-linear heat conduction
% problem. This solution will be used as a test example for the inverse 
% solver.
%
% We monitor the convergence. A couple of iterations is enough to get
% good accuracy. Here we do more iterations to get the convergence graph.
%
%======================================================================%

 tic,
 IterationError=zeros(1,0);T=zeros(size(X));
 for i=1:15,
   [kappa,Ap]=ThermalModel(X,Z,T);
   [Tnew,Q0]=DirectThermalSolve( x , z , kappa , Ap , T0 , Qm );
   IterationError(i)=norm(Tnew-T,'fro');
   T=Tnew;
 end
 Te=Tnew;
 semilogy(IterationError,'LineWidth',1.5);
 ylabel('Step size: ||T^{(k)}-T^{(k-1)}||_F','FontSize',14);
 xlabel('Iteration number: k','FontSize',14)
 toc
 
 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'The convergence of the fixed point iterations.\n\n');
 input('Press return to continue');
 
% Use the computed solution to display the Thermal conductivity and 
% Heat production in the Litosphere

 close all
 [kappa,Ap]=ThermalModel(X,Z,T,'on');
 
 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'The Thermal conductivity and heat production.\n\n');
 input('Press return to continue');

% Display the temperature distribution T(x,z) for the numerical example
% created above. This is the solution of the direct problem.

 close all
 mesh(X(1:5:end,1:5:end)/10^3,Z(1:5:end,1:5:end)/10^3,T(1:5:end,1:5:end))
 xlabel('Horizontal: x [ km ]','FontSize',14);
 ylabel('Depth: z [ km ]','FontSize',14);
 zlabel('Temperature: T(x,z) [ ^oC ]','FontSize',14)
 
 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'The Temperature distribution T(x,z) for the numerical example.\n\n');
 input('Press return to continue');
 
 

%=========================================================================%
%
% By running the above Matlab lines we have effectively created a simulated
% problem with a known solution. We have both T(x,0)=T0 and kT_x(x,0)=Q0
% for a known solution Te(x,z). Now we use this known solution to test the 
% solution method for the inverse problem. We add simulated noise to the 
% data. 
%
% Tikhonov regularisation is used to stabilize the computations. The
% algorithm implements Tikhonov regularization by an iterative procedure
% that uses the linear system of equations produced by the finite
% difference discretization of the direct problem.
%
%Also produce an L-curve and a plot of Error 
% vs Lambda. The problem is solved by computing the matrix representation
% of the operator explicitly so can also compute the singular values. 
%
%=========================================================================%


%
% Add errors to both T0, and Q0. This should ideally be of the same
% order of magnitude as the noise in a real measurement. If changing here
% you also need to change the regularization parameter later.
%
 T0err=T0+randn(size(T0))*5e-3;
 Q0err=Q0+randn(size(Q0))*5e-5;
 

% First illustrate the noisy Cauchy data. 

 close all
 plot(x/10^3,-10^3*Q0,'b--',x/10^3,-10^3*Q0err,'k-','LineWidth',1.5);
 xlabel('Horizontal Coordinate: x [ km ]','FontSize',14);
 ylabel('Heat-flux: Q_0 [ mW/m ]','FontSize',14)
 axis([0 400 65 73])
 figure
 plot(x/10^3,T0,'b--',x/10^3,T0err,'k-','LineWidth',1.5);
 xlabel('Horizontal Coordinate: x [ km ]','FontSize',14);
 ylabel('Surface temperature: T_0 [ ^oC ]','FontSize',14)
 axis([0 400 9.8 10.2])
 
  
 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'The noisy simulated Cauchy data T0 and Q0.\n\n');
 input('Press return to continue');
 


% Solve the problem using a specific Lambda. Since the problem is
% non-linear we have to again use fixed point iteration. The initial
% guess is zero temperature in the domain. This is bad and causes the 
% initial error to be huge. Still only a couple of iterations is needed.
% In each fixed point iteration we use the previous solution as initial
% guess to the PCG iterations. 
%
% Increase Lambda for more regularization and decrease for less. 
%

 Lambda=6.3e-3;
 T=zeros(size(Te));Qtik=zeros(size(Q0err));
 IterationError=zeros(1,0);
 tic
 for i=1:10,
  [K,Ap]=ThermalModel(X,Z,T);
  [Qtik]=LinearTikhonovSolve(x,z,K,Ap,T0err,Q0err,Lambda,Qtik);
  [Tnew]=DirectThermalSolve(x,z,K,Ap,T0err,Qtik);
  IterationError(i)=norm(T-Tnew,'fro');
  fprintf(1,'Iteration:  %i  Error: %e  ',i,IterationError(i));
  T=Tnew;toc
 end; 
 


 close all
 semilogy(IterationError,'LineWidth',1.5);
 ylabel('Step size: ||T^{(k)}-T^{(k-1)}||_F','FontSize',14);
 xlabel('Iteration number: k','FontSize',14)
 
 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'The convergence for the fixed point iterations for the inverse problem.\n\n');
 input('Press return to continue');
 

 plot(x/10^3,10^3*Qm,'b--',x/10^3,10^3*Qtik,'k-','LineWidth',1.4)
 xlabel('Horizontal Coordinate: x [ km ]','FontSize',14);
 ylabel('Heat-flux: Q_m and Q_{tik} [ mW/m ]','FontSize',14)
 axis([0 400 20 40])
 
 fprintf(1,'\n---------------------------------------------------\n');
 fprintf(1,'The Tikhonov solution Qtik and also the exact solution Qm.\n\n');
 input('Press return to continue');
   
 
