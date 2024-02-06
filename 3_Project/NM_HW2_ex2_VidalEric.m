%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window
%source('mystartdefaults.m'); % Contains SI physical constants
mystartdefaults; % Contains SI physical sonstants

tic

%% UNITS FOR FREE ELECTRONS IN VACUUM
recipunit = 1.0E+10; % Å^-1, unit of the reciprocal length
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm)) / qel;

%% COMPUTING THE 1D MODEL OF AN ELECTRON SCATTERED BY A QUANTUM JUNCTION
% THIS IS BASED IN SOLVING THE TIME-INDEPENDENT SCHRODINGER EQUATION
% BY THE LOCALIZED GREEN'S FUNCTION METHOD FOR A PARTICLE IN A POTENTIAL WELL

% The quantum junction is two potential barriers of height U(x')=0.2eV for
% x' in [0,15]Å and [65,80]Å, and U(x)=0eV elsewhere
% The whole x range is [-20,100]Å
% Then an electric field is applied, so it modifies the potential,
% resulting in a biased potential.

% This resembles the order of magnitude of a heterostructure analog to
% resonant tunneling devices in GaAs-AlGaAs

%% GRID DEFINITION
% Discretized grid for the potential U(x)
x_step = 0.5    % Å, Step for the grid 
xmin = 0        % Å, Minimum x value
xmax = 80       % Å, Maximum x value
x_U = (xmin + x_step/2):x_step:(xmax - x_step/2); % Å, Discretized grid for x
% Remaining grid
x_neg = (-20 + x_step/2):x_step:(xmin - x_step/2); % Å, Discretized grid for x
x_pos = (xmax + x_step/2):x_step:(100 - x_step/2); % Å, Discretized grid for x
% Total grid
xx = [x_neg, x_U, x_pos]; % Å, Discretized grid for x

% Potential U(x) for the quantum junction
Ux = BarrierPotential(x_U, 0, 15, 0.2) + BarrierPotential(x_U, 65, 80, 0.2);
Ufull = [zeros(1,length(x_neg)), Ux, zeros(1,length(x_pos))];

% Energy discretization
E_step = 0.0005  % eV, Step for the energy
Emin = 0.0       % eV, Minimum energy
Emax = 0.3       % eV, Maximum energy
EE = (Emin + E_step/2):E_step:(Emax - E_step/2); % eV, Discretized energy
% E length
numE= length(EE);
% Recombination time (lifetime) parameter and damping factor
recombT = 1.0E-9; % ns
damping=(hbar * 2*pi / recombT) / qel;

% Potential bias due to an electric field aplied to the heterojunction for x>0
% Liniar potential with a slope of -0.1eV/Å and width of the heterostructure
Ebias = -0.1; % eV/Å, Slope of the potential bias
for i=1:length(x_U)
    Ubiased(i) = Ux(i) + Ebias*x_U(i)/(xmax-xmin);
end
Ufull_biased = [zeros(1,length(x_neg)), Ubiased, Ebias*ones(1,length(x_pos))];

% Plot the biased potential U(x)
figure(1)
plot(xx, Ufull, 'b', xx, Ufull_biased, 'r')
% title('Potential U(x) for the quantum junction')
xlabel('x [Å]')
ylabel('U(x) [eV]')
ylim([-0.15,0.35])
fontsize(gca, 22,'points')   % 'pixels', 'centimeters', 'inches'
set(gca,'Box','on');
set(gca,'linewidth',1);
legend('U(x)', 'U(x) + U_{biased}(x)')
grid on
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);


%% PERTURBATION
% According to equation (4.13) for solving the step perturbation with
% V' at the end is what we had in the total potential - 0th order solution:
Uperturbation = Ubiased - Ebias*ones(1,length(x_U));
Ucheck = [zeros(1,length(x_neg)), Uperturbation, zeros(1,length(x_pos))];
% where the 0th order solution is the step potential

% Notice that from -\infty until 0 it is 0, and from 80 until \infty it is 0
% so (eq. A.65) can be applied which come next

%% SANITY CHECK
% Plot it
figure(2)
plot(xx, Ucheck, 'b')
% title('Perturbation U(x) for the quantum junction')
xlabel('x [Å]')
ylabel('U(x) [eV]')
ylim([0,0.35])
fontsize(gca, 22,'points')   % 'pixels', 'centimeters', 'inches'
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);

%% WAVE FUNCTIONS
% It comes from the dispersion relations binding E to the wavevector
% is not conserved across the step. Leading to eq. 58 in the book of NM
commonTerm = 2*elm*qel/hbar^2;
% for x<0
k0min = real(sqrt((Emin + 1i*damping) * commonTerm));
wavelenght0max = 2*pi*(1E+10)/k0min % Å, Shortest Wavelength
k0max = real(sqrt((Emax + 1i*damping) * commonTerm));
wavelenght0min = 2*pi*(1E+10)/k0max % Å, Longest Wavelength
% for x>0
k1min = real(sqrt((Emin - Ebias + 1i*damping) * commonTerm));
wavelenght1max = 2*pi*(1E+10)/k1min % Å, Shortest Wavelength
k1max = real(sqrt((Emax - Ebias + 1i*damping) * commonTerm));
wavelenght1min = 2*pi*(1E+10)/k1max % Å, Longest Wavelength


%% LOCALIZED GREEN'S FUNCTION METHOD
% The Green's function for the Schrodinger equation is given by
% (eq. A.76) and (eq. A.77) in the book of NM

% Production of R(E), T(E) and A(E) curves for the biased quantum junction
% Same conditions as in Exercise 1, except step which is the bias now.
[RR, TT, AA] = RTA(Ebias, EE, damping, x_U, Uperturbation, x_step, ekinscale);


%% ELECTRON AT THE FERMI LEVEL
% Set E0=0.01eV and explore the effect of tuning the bias
% U1 in [-0.2, 0.2]eV
U0min = -0.2; % eV
U0max = 0.2; % eV

for i=1:length(x_U)
    Ubiased(i) = Ux(i) + U0min*x_U(i)/(xmax-xmin);
end
U0min_biased = [zeros(1,length(x_neg)), Ubiased, U0min*ones(1,length(x_pos))];

for i=1:length(x_U)
    Ubiased(i) = Ux(i) + U0max*x_U(i)/(xmax-xmin);
end
U0max_biased = [zeros(1,length(x_neg)), Ubiased, U0max*ones(1,length(x_pos))];
% Plot the biased potential U(x)
figure(4)
plot(xx, Ufull, 'b', xx, U0min_biased, 'r', xx, U0max_biased, 'g')
% title('Potential U(x) for the quantum junction')
xlabel('x [Å]')
ylabel('U(x) [eV]')
ylim([-0.25,0.45])
fontsize(gca, 22,'points')   % 'pixels', 'centimeters', 'inches'
set(gca,'Box','on');
set(gca,'linewidth',1);
legend('U(x)', 'U(x) + U_{biased}(x) [U_{1}(x)=-0.2eV]', 'U(x) + U_{biased}(x) [U_{1}(x)=0.2eV]', 'Location', 'Best')
grid on
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);


%% RTA FOR THE ELECTRON AT THE FERMI LEVEL
% Fine U1 step of 0.0005eV, compute current-bias curve in this range of bias
% Electrical resistance, diode? Resonant-tunneling diode!
setE=0.01; % eV, E0

EE = (Emin + E_step/2 + U0min):E_step:(Emax - E_step/2 + U0min); % eV, Discretized energy

for i=1:length(EE)
    for j=1:length(x_U)
        perturbation(j) = EE(i)*x_U(j)/(xmax-xmin) - EE(i);
    end
    Upertu=Ux+perturbation;
    [RR(i), TT(i), AA(i)] = RTA_iter(EE(i), setE, damping, x_U, Upertu, x_step, ekinscale);
end


% Plot of R(E), T(E) and A(E)
% Note that the x axis is reversed
figure;
plot(-EE, RR, 'LineWidth', 2);
hold on;
plot(-EE, TT, 'LineWidth', 2);
plot(-EE, AA, 'LineWidth', 2);
% title('R(E), T(E) and A(E) for the Biased Quantum Junction', 'fontsize', 26);
xlabel('E (eV)','FontSize',26);
ylabel('R(E), T(E), A(E)','FontSize',26);
fontsize(gca, 22,'points')   % 'pixels', 'centimeters', 'inches'
ylim([0,1]);
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on;
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);
legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');


% Plot of R(E), T(E) and A(E)
% Note that the x axis is reversed
figure;
plot(-EE, RR, 'LineWidth', 2);
hold on;
plot(-EE, TT, 'LineWidth', 2);
plot(-EE, AA, 'LineWidth', 2);
% title('R(E), T(E) and A(E) for the Biased Quantum Junction', 'fontsize', 26);
xlabel('E (eV)','FontSize',26);
ylabel('R(E), T(E), A(E)','FontSize',26);
fontsize(gca, 22,'points')   % 'pixels', 'centimeters', 'inches'
xlim([0.06,0.07]);
ylim([0,0.2]);
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on;
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);
legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');


toc