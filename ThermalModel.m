% This is a simple model of the layers in the ground and their
% thermal properties. Can also display results in convinient ways.
%
%
%
%
function [K,A]=ThermalModel(X,Z,T,DisplayGraph)


%
% Extract vectors x and z from the matrices X and Z.
%
 x=X(1,:)';if x(1)==x(2),x=X(:,1);,end
 z=Z(1,:)';if z(1)==z(2),z=Z(:,1);,end

 [LayerDepth]=SetLithosphericStructure( x );
 [K]=SetThermalConductivity( X , Z , T , LayerDepth  );
 [A]=SetHeatProduction( X , Z , LayerDepth  );

  %
  % If the DisplayGraphs command is given we produce plots for the paper.
  % 
   if nargin==4,% Display thermal conductivity and temperature
      figure
      DisplayLithosphericStructure( x , z , LayerDepth , K , 'Thermal Conductivity: k [W/m ^oC]' );
      figure
      DisplayLithosphericStructure( x , z , LayerDepth , A , 'Heat Production: A_p [W/m^3]' );
    end

end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 function [K]=SetThermalConductivity( X , Z , T , LayerDepth  )
 
 [M ,N]=size(X);      % Number of gridpoints in z- and x-direction
  K=zeros(M,N);
  %
  % Default is values for the mantle.
  %
   K=2.5*ones(M,N)./(1-2.5e-4*T);
  %
  % Lower crust layer
  %
   ind=find(Z<ones(M,1)*LayerDepth(:,3)');
   K(ind)=2.6./(1+1.0e-4*T(ind));
 
  %
  % Upper crust layer
  %
    ind=find(Z<ones(M,1)*LayerDepth(:,2)');
    K(ind)=3.0./(1+1.5e-3*T(ind)); 
  %
  % Sedimemt layer
  %
    ind=find(Z<ones(M,1)*LayerDepth(:,1)');
    K(ind)=3.0./(1+0*T(ind)); 

 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [A]=SetHeatProduction( X , Z , LayerDepth )
 
 [M ,N]=size(X);      % Number of gridpoints in z- and x-direction
  
  %
  % Set to default equal to lowest layer. Then overwrite. 
  %
   A=0.02e-6*ones(M,N);
  %
  % Lower crust layer
  %
   ind=find(Z<ones(M,1)*LayerDepth(:,3)'); 
   A(ind)=0.45e-6;                
  %
  % Upper crust layer
  %
   ind=find(Z<ones(M,1)*LayerDepth(:,2)'); 
   A(ind)=2.0e-6;
   
  %
  % Sediment layer
  % 
   ind=find(Z<ones(M,1)*LayerDepth(:,1)'); % sediment
   A(ind)=1.50e-6;                

  
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [LayerDepth]=SetLithosphericStructure( x )

    
  %
  % Set the depth of the different litospheric layers as a function of 
  % x. The length of the physical model is L1.
  % 
     
   LayerDepth=zeros(length(x),4);
   LayerDepth(:,1)=(5-2*cos(2*pi*x/max(x)))*1e3;    % sediment depth     [m]
   LayerDepth(:,2)=(15+3*cos(2*pi*x/max(x)))*1e3;   % upper crust depth  [m]
   LayerDepth(:,3)=(35+5*cos(2*pi*x/max(x)))*1e3;   % moho depth         [m]
   LayerDepth(:,4)=(120+10*cos(2*pi*x/max(x)))*1e3; % lithosphere depth  [m]

 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DisplayLithosphericStructure( x , z , LayerDepth , K , TitleText )

  % Reverse the z-coordinate so its 0 to -L1. Also change to km.
  z=-z/1e3;x=x/1e3;LayerDepth=-LayerDepth/1e3; 
  clf,
  colormap(jet);
  pcolor(x(1:2:end),z(1:2:end),K(1:2:end,1:2:end));shading interp;hold on
  colorbar;
  title(TitleText,'FontSize',14)
  xlabel('Horizontal Position x [km]','FontSize',14)
  ylabel('Depth [km]','FontSize',14)
  set(get(gcf,'CurrentAxes'),'YDir','normal');
  plot([min(x),max(x),max(x),min(x),min(x)],[min(z),min(z),max(z),max(z),min(z)],'k','LineWidth',1.5);
  borderLength=0.05*(max(z)-min(z));
  axis([min(x)-borderLength max(x)+borderLength min(z)-borderLength max(z)+borderLength]);
  hold on
  plot(x,LayerDepth(:,1),'k','LineWidth',1.4);
  plot(x,LayerDepth(:,2),'k','LineWidth',1.4);
  plot(x,LayerDepth(:,3),'k','LineWidth',1.4);
  plot(x,LayerDepth(:,4),'k','LineWidth',1.4);
end