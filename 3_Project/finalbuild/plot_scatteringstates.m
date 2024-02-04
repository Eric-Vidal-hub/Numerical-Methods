%% Script for plotting scattering states of selected energies
clearvars
Energies=[0.01, 0.04325, 0.09525, 0.16325, 0.24075, 0.06]; 
%% Plot flags
plotpotential = 'n';
plotpeturbsoln = 'n';
plotfullwavefn = 'n';

%% units and constants
mystartdefaults
recipunit = 1e10;          % Reciprocal space unit = 1e10 1/m
% Constant to convert k (A^2) to Energy (eV)
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;
fprintf('\n ekinscale = %d \n',ekinscale);

%% Discretisation
spatial_disc

%% define peturbation
UL = 0;       % eV Potential of reference system if x<0
UR = -0.0;    % eV Potential of reference system if ><0

U1 = 0.2;       % eV                   % Height of 1st barrier
U2 = 0.2;       % eV                   % Height of 2nd barrier
peturb_init

lifetime = 1e-9;
gamma = hbar*(2*pi/lifetime)/qel;

%% Calculate states
states = zeros(L1,length(Energies));
legendentries{1} = 'Potential * ekinscale';
for ii=1:length(Energies)
    E0 = Energies(ii);             % eV
    calc_peturbsoln
    calc_propogation
    states(:,ii) = abs(phis).^2;
    legendentries{ii+1} = sprintf('E0 = %2.3d',E0);
end
%% Plot states
stateplot = figure;
plot(x,Vfull*ekinscale,'black')
hold on
ind = sort(1:length(Energies)-1,'descend');
for ii=ind
    plot(x,states(:,ii) - 0.1*(5-ii) - 0.5,'blue','LineWidth',1.5)
end
plot(x,states(:,end)-0.75,'red','LineWidth',1.5)

title('Resonant Modes','FontSize',15)
grid on
legend(legendentries,'FontSize',15,'location','eastoutside')