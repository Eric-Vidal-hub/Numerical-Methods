% Script to define and solve wave within peturbation and 
% propogate solution to full grid outside peturbation, 
% with peturbation defined by calling initialconds.m
%% Defining incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 1;                % Normalisation factor

k0 = sqrt((E0-UL+(1i*gamma))/ ekinscale); % Eq. A.58
k1 = sqrt((E0-UR+(1i*gamma))/ ekinscale); % Eq. A.58

rb = (k0-k1)/(k1+k0); % Eq. A.64
tb = (2*k0)/(k1+k0); % Eq. A.64

phi0=tb*exp(1i*k1*Ux);% Eq. A.65
%% Solve system 

% Greens Function 
G1 = Greenfunc(Ux,Ux,k0,k1);      % all Ux>=0 =>  k1,k1

% Scattering Matrix - Eq 4.51
T = eye(L2,L2) - xres*G1*diag(W);   % Scattering Matrix
phip = T\(phi0);           % Solve System

%% plot
if(plotpeturbsoln=='y')
fig1 = figure;
hold on
grid on
plot(Ux, abs(phip).^2)
plot(Ux, W)
title(sprintf('Solution within peturbation, E0=%5.3d',E0),'FontSize',15)
legend('Wavefunction','Potential','FontSize',15)
end