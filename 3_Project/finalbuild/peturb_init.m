%% Defines the peturbation for a pre-defined spatial grid
% calls upon naming conventions from spatial_disc.m

%% reference system

Vreffull = zeros(length(x),1);
Vref = zeros(length(Ux),1);

for i = 1:L2
  if Ux(i) < 0
    Vref(i) = UL;
  end
  if Ux(i) >= 0
    Vref(i) = UR;
  end
end

for i = 1:L1
  if x(i) < 0
    Vreffull(i) = UL;
  end
  if x(i) >= 0
    Vreffull(i) = UR;
  end
end

%% Defining perturbation 
V = zeros(L2,1);
modE = -(UR-UL)/(Uxmax2-Uxmin1);

for i = 1:L2
  if Ux(i) >= Uxmin1 && Ux(i) <= Uxmax1
    V(i) = U1;  
  end
  if Ux(i) >= Uxmin2 && Ux(i) <= Uxmax2
    V(i) = U2;
  end
end
%% applying bias
V = V - modE*Ux; 
Vfull = zeros(L1,1);
for i = 1:L1
  if x(i) >= Uxmin1 && x(i) <= Uxmax1
    Vfull(i) = U1;
  end
  if x(i) >= Uxmin2 && x(i) <= Uxmax2
    Vfull(i) = U2;
  end
end

%% Peturbation with respect to reference
for i=1:L1
    if x(i)>=0 && x(i)<Uxmax2
        Vfull(i) = Vfull(i) - modE*x(i);
    end
    if x(i) >= Uxmax2
        Vfull(i) = UR;
    end
end

W = V - Vref;
Wfull = Vfull - Vreffull;

%% Plot
if(plotpotential=='y')
    potplot=figure;  
    plot(x,Vreffull,'LineWidth',1.5); hold on;
    plot(x,Vfull,'LineWidth',1.5); hold on;
    plot(x,Wfull,'LineWidth',1.5); hold on;
    
    % xlim([xmin,xmax]);
    % ylim = [0, 03];

    xlabel('x (angstrom)','FontSize',15);
    ylabel('(eV)','FontSize',15);
    axis([xmin, xmax, -0.15, 0.5])

    legendU{1} = 'Vreffull';   
    legendU{2} = 'Vfull';
    legendU{3} = 'Wfull';
    lgd=legend(legendU,'location','eastoutside');  
  
    grid on;
end

%% Energy limits due to sampling grid
lambdamax = 2*(max(Ux) - min(Ux));
lambdamin = xres;

Emin = ((hbar^2)/(2*elm*(lambdamax*directunit)^2)) + UL*qel;
Emax = ((hbar^2)/(2*elm*(lambdamin*directunit)^2)) + UL*qel + UR*qel;

fprintf('\n Spatial discretisation imposes limits on incident electron \n')
fprintf('\n      Minimum energy E_a = %5.3d eV\n      Maximum energy E_b = %5.3d eV\n',Emin/qel, ...
    Emax/qel)