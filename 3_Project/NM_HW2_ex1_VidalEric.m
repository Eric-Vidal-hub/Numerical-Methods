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
E_step = 0.0005  % eV, Step for the energy
Emin = 0.0       % eV, Minimum energy
Emax = 0.3       % eV, Maximum energy
EE = (Emin + E_step/2):E_step:(Emax - E_step/2); % eV, Discretized energy
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

% Plot of R(E), T(E) and A(E)
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

% Plot only of A(E)
figure;
plot(EE, AA, 'LineWidth', 2);
title('A(E) for the Quantum Junction', 'fontsize', 26);
xlabel('E (eV)','FontSize',18);
ylabel('A(E)','FontSize',18);
ylim([0,0.03]);
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on;
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);

%% RESONANT TUNNELING DIODE (RTD) MODEL
% The RTD model is a quantum junction with a resonant level in the middle
% of the potential barriers

% Determining the resonant Energies from the past Figure
resonanceE = [0.01075, 0.04325, 0.09525, 0.16325, 0.23975]
nonresonantE = (resonanceE(3) + resonanceE(4)) / 2
% Concatenate the resonant and nonresonant energies
resoE = [resonanceE, nonresonantE]

[~, nreso] = size(resoE);
wavefunctions = zeros(nreso, length(x_neg) + length(x_U) + length(x_pos));
Wmat = zeros(length(x_U), length(x_U));
for ii = 1:nreso
    k1 = sqrt((resoE(ii) + 1i * damping) / ekinscale);
    for jj = 1:length(x_U)
        phiok(jj) = exp(1i * k1 * x_U(jj));
        Wmat(jj, jj) = x_step * Ux(jj) / ekinscale;
        for kk = 1:length(x_U)
            Go(jj, kk) = GreensFun(0, x_U(jj), x_U(kk), resoE(ii), damping, ekinscale);
        end
    end
    scatt = eye(length(x_U)) - Go * Wmat;
    phisol = scatt \ (phiok.'); 
    for jj = 1:(length(x_neg) + length(x_U) + length(x_pos))
        outside_sol_val = ExtraTerm(xx(jj), x_U, Ux, x_step, phisol, 0, resoE(ii), damping, ekinscale);
        if ii < nreso
            wavefunctions(ii, jj) = exp((xx(jj)) * 1i * k1) + outside_sol_val;
        else
            nonresophi(jj) = exp((xx(jj)) * 1i * k1) + outside_sol_val;
        end
    end
end
resoprobas = abs(wavefunctions).^2;
nonresoproba = abs(nonresophi).^2;

% Plot of the scaled resonant wavefunctions for readability and comparison
figure;
hold on;
plot(xx, (100/0.2)*Ufull, 'LineWidth', 1, 'Color', 'black');
plot(xx, (100/4.)*resoprobas(1,:), 'LineWidth', 1);
plot(xx, (100/27.)*resoprobas(2,:), 'LineWidth', 1);
plot(xx, (100/112.)*resoprobas(3,:), 'LineWidth', 1);
plot(xx, (100/30.)*resoprobas(4,:), 'LineWidth', 1);
plot(xx, (100/6.)*resoprobas(5,:), 'LineWidth', 1);
% How to set color gray
plot(xx, (100/4.)* nonresoproba, 'LineWidth', 1, 'Color', [0 0 0]+0.5);
title('Scaled Resonant Wavefunctions for the Quantum Junction', 'fontsize', 26);
xlabel('x (Å)','FontSize',18);
ylabel('|\psi(x)|^2','FontSize',18);
legend('U(x)', '|\psi_1(x)|^2', '|\psi_2(x)|^2', '|\psi_3(x)|^2', '|\psi_4(x)|^2', '|\psi_5(x)|^2', '|\psi_6(x)|^2', 'Location', 'Best');
ylim([0,125]);
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on;
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);


%% CRUDE APPROACH TO AN APPLIED BIAS
% Change the potential U(x) for the quantum junction by lowering the
% potential barrier at the right side of the junction to 0.1eV
Ux = BarrierPotential(x_U, 0, 15, 0.2) + BarrierPotential(x_U, 65, 80, 0.1);
% Production of R(E), T(E) and A(E) curves for the quantum junction
energyStep = 0;      % constant background
for i = 1:numE
    % For each energy level, calculate the reflection (R), transmission (T), and absorption (A) coefficients
    % using the RTA function with a parallelized for loop
    [RR(i), TT(i), AA(i)] = RTA(energyStep, EE(i), damping, x_U, Ux, x_step, ekinscale);
end

% Plot of R(E), T(E) and A(E)
figure;
plot(EE, RR, 'LineWidth', 2);
hold on;
plot(EE, TT, 'LineWidth', 2);
plot(EE, AA, 'LineWidth', 2);
title('R(E), T(E) and A(E) for the Modified Quantum Junction', 'fontsize', 26);
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

% Plot only of A(E)
figure;
plot(EE, AA, 'LineWidth', 2);
title('A(E) for the Modified Quantum Junction', 'fontsize', 26);
xlabel('E (eV)','FontSize',18);
ylabel('A(E)','FontSize',18);
ylim([0,0.005]);
set(gca,'Box','on');
set(gca,'linewidth',1);
grid on;
% Set the color of the axes and the grid to a light gray
set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
% Set the background color of the figure to a darker gray
set(gcf, 'Color', [0.7 0.7 0.7]);
% set(gca,'TickLength',[0.03, 0.02]);

toc