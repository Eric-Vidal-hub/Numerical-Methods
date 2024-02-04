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
% x' in [0,15]Å and [65,80]Å, and U(x)=0eV elsewhere
% The whole x range is [-20,100]Å

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
title('Potential U(x) for the quantum junction')
xlabel('x [Å]')
ylabel('U(x) [eV]')
set(gca,'Box','on');
set(gca,'linewidth',1);
legend('U(x)', 'U(x) + U_{bias}(x)')
grid on
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);


%% LOCALIZED GREEN'S FUNCTION METHOD
% The Green's function for the Schrodinger equation is given by
% (eq. A.76) and (eq. A.77) in the book of NM

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
title('Perturbation U(x) for the quantum junction')
xlabel('x [Å]')
ylabel('U(x) [eV]')
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);



toc