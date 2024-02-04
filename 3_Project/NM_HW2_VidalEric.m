%% START DEFAULT COMMANDS OF *.m SCRIPTS
close all;       % Close all figures if any
clear all;       % Clear all variables/functions in memory
clc;             % Clear screen in the command window
%source('mystartdefaults.m'); % Contains SI physical sonstants
mystartdefaults; % Contains SI physical sonstants

tic

%% UNITS FOR FREE ELECTRONS IN VACUUM
recipunit = 1.0E+10; % Å^-1, unit of the reciprocal length
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm)) / qel;

%% COMPUTING THE 1D MODEL OF AN ELECTRON SCATTERED BY A QUANTUM JUNCTION
% THIS IS BASED IN SOLVING THE TIME-INDEPENDENT SCHRODINGER EQUATION
% BY THE LOCALIZED GREEN'S FUNCTION METHOD FOR A PARTICLE IN A POTENTIAL WELL

% The quantum junction is two potential barriers of height U(x')=0.2eV for
% x' in [0,15]Å and [65,80]Å, and U(x)=0eV elswhere
% The whole x range is [-20,100]Å

% This resembles the order of magnitude of a heterostructure analog to
% resonant tunneling devices in GaAs-AlGaAs

%% GRID DEFINITION
% Discretized grid for the potential U(x)
x_step = 0.5;   % Å, Step for the grid 
xmin = 0;       % Å, Minimum x value
xmax = 80;      % Å, Maximum x value
x_U = (xmin + x_step/2):x_step:(xmax - x_step/2); % Å, Discretized grid for x
% Remaining grid
x_neg = (-20 + x_step/2):x_step:(xmin - x_step/2); % Å, Discretized grid for x
x_pos = (xmax + x_step/2):x_step:(100 - x_step/2); % Å, Discretized grid for x
% Total grid
xx = [x_neg, x_U, x_pos]; % Å, Discretized grid for x


% Potential U(x) for the quantum junction
Ux = BarrierPotential(x_U, 0, 15, 0.2) + BarrierPotential(x_U, 65, 80, 0.2);
Ufull = [zeros(1,length(x_neg)), Ux, zeros(1,length(x_pos))];
% Create a new figure
figure;
plot(xx, Ufull, 'LineWidth', 2);
title('Barrier Potential', 'fontsize', 26);
xlabel('x (Å)','FontSize',18);
ylabel('U (eV)','FontSize',18);
ylim([0,0.22]);
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on;
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);

%% LOCALIZED GREEN'S FUNCTION METHOD
% Energy discretization
E_step = 0.0005; % eV, Step for the energy
Emin = 0.0;      % eV, Minimum energy
Emax = 0.3;      % eV, Maximum energy
EE = (Emin + E_step/2):E_step:(Emax - E_step/2) % eV, Discretized energy
% E length
numE= length(EE);
% Recombination time (lifetime) parameter and damping factor
recombT = 1.0E-9; % ns
damping=(hbar * 2*pi / recombT) / qel;

% Production of R(E), T(E) and A(E) curves for the quantum junction
energyStep = 0;      % constant background
for i = 1:numE
    % For each energy level, calculate the reflection (R), transmission (T), and absorption (A) coefficients
    % using the RTA function with a parallelized for loop
    [RR(i), TT(i), AA(i)] = RTA(energyStep, EE(i), damping, x_U, Ux, x_step, ekinscale);
end

% Create a new figure
figure;
plot(EE, RR, 'LineWidth', 2);
hold on;
plot(EE, TT, 'LineWidth', 2);
plot(EE, AA, 'LineWidth', 2);
title('R(E), T(E) and A(E) for the Quantum Junction', 'fontsize', 26);
xlabel('E (eV)','FontSize',18);
ylabel('R(E), T(E), A(E)','FontSize',18);
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

toc


