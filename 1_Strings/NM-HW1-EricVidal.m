% Task: NM-HW1 2.9.5 Damped Vibrating String
% Deadline: October 4th, 2023 at 12:00.
% Author: Eric Vidal Marcos
% Group: 2 (QuanTEEM)
% Email: ericvidalmarcos@gmail.com
% Supervisors: Dr. Alain Dereux & Dr. Matthieu Sala
% Octave Version: 8.3.0

% Clear workspace
clear all   % Clear all variables
close all   % Close all variables
clc         % Clear command window

% Plot commands
graph_pinch = 'yes';            % Plot initial conditions of f and g
graph_eigenfunctions = 'yes';   % Plot the eigenfunctions
graph_f_approx = 'yes';         % Plot the comparison of f and approx f
graph_time = 'yes';             % Plot the evolution of the wave
graph_movie = 'yes';            % Plot the movie of the wave evolution

% Pinch commands
minpinch = 1;                   % Minimum number of pinches
maxpinch = 2;                   % Maximum number of pinches
% To set just one pinch, minpinch = maxpinch

% 1. Input data and discretization parameters
L = 0.328           % String length [m]
mu = 0.660E-3       % String mass per unit length (mu = M / L) [kg/m]
F = 55.0            % Tension force [N]
cc = sqrt(F/mu)     % Speed of wave propagation [m/s] (eq. 2.48)
nmax = 20           % Maximum number of modes (to be adjusted for the good convergence of the approximation)
npt = 201           % Number of sampling points used in the approximation
nst = 13            % Number of sample times
% movie=200           % Number of frames for movie development phase

gg = 100.0                  % Gamma coefficient [Hz]
dampf = 2 * gg / (cc ^ 2)   % Damping factor [s/m] (eq. 2.86)
tmax = (2 / gg) * log(10)   % Maximum time [s] (my pdf eq. 3)
tmin = 0;                   % Minimum time [s]
tstep = (tmax - tmin) / (nst - 1);  % Time step [s]
for i=1:nst
  tt(i) = tmin + (i - 1) * tstep;    % Discretized time [s]
end

% 2. Array of discretized value of x
xmin = 0;           % Minimum abscissa  
xmax = L;           % Maximum abscissa  
xstep = (xmax - xmin) / (npt - 1);  % Discretization step [m]
for i = 1:npt
  x(i) = xmin + (i - 1) * xstep;    % Discretized abscissas  
end

% 3. Compute numerical arrays for both pinches and to be used as initial conditions
for pinch = minpinch:maxpinch
    if (pinch == 1);  % eq. (2.44)
        p=L/2;   % Pinched point
        z = 1;   % Amplitude of the wave
        for (i=1:npt)   % Initial deformation for Pinch 1
            if (x(i) <= p)
                f(pinch, i) = x(i) * z / p;
            else
                f(pinch, i) = z * (x(i) - L) / (p-L);
            end
        end
        for i=1:npt     % Initial distribution of velocity field
            g(pinch, i) = 0;
        end
    elseif (pinch == 2);  % eq. (2.45)
        p1 = L/4.;    % Pinched point 1
        z1 = 1;       % Amplitude of the wave 1

        p2 = (3./4.) * L; % Pinched point 2
        z2 = -z1;       % Amplitude of the wave 2
        for (i=1:npt)   % Initial deformation for Pinch 2
            if (xmin <= x(i)) && (x(i) <= p1)
                f(pinch, i) = x(i) * z1 / p1;
            elseif (p1 < x(i)) && (x(i) <= p2)
                f(pinch, i) = (z1 - z2) / (p1 - p2) * (x(i) - p1) +z1;
            elseif (p2 < x(i)) && (x(i) <= xmax)
                f(pinch, i) = z2 * (x(i) - L) / (p2 - L);
            end
        end
        for i=1:npt     % Initial distribution of velocity field
            g(pinch, i) = 0;
        end
    end
end

% Plot the initial conditions of f and g for each pinch
if (strcmp(graph_pinch, 'yes') == 1);
    for pinch = minpinch:maxpinch;
        p1 = plot(x/L, f(pinch, :), ";f(x) [m];", 'LineWidth', 2);
        hold on
        p1 = plot(x/L, g(pinch, :), ";g(x) [m/s];", 'LineWidth', 2);
        hold off
        xlim([-(xmax / L) * 0.1, (xmax / L) * 1.1]);
        ylim([min(f(pinch,:)) - 0.1, max(f(pinch,:)) + 0.1]);
        xlabel('x/L', 'FontSize', 14);
        ylabel('y(x)', 'FontSize', 14);
        title(strjoin({'Pinch', num2str(pinch)}, ' '), 'FontSize', 14);
        legend("location", "northeast", 'FontSize', 14);
        waitfor(p1)   % Stops the code execution until you close the plot
    end
end

% 4. Compute the eigenfunctions and plot
for n = 1:nmax
    for i = 1:npt
        phi(n, i) = sqrt(2 / L) * sin(n * pi * x(i) / L);   % Eigenfunctions (eq. 2.75)
    end
end

% GRAPH OF EIGENFUNCTIONS
if (strcmp(graph_eigenfunctions, 'yes') == 1);
    for n = 1:nmax
        p2 = plot(x/L, phi(n, :), 'LineWidth', 1);      % For eigenfunctions labeling include strjoin({';Mode ', num2str(n), ';'}, ''),
        hold on;
    end
    xlim([-(xmax / L) * 0.1, (xmax / L) * 1.1]);
    ylim([-sqrt(2/L) * 1.1,sqrt(2/L) * 1.1]);
    xlabel('x/L', 'FontSize', 14);
    ylabel('{\Phi(x)} [m]', 'FontSize', 14);
    title('Eigenfunctions', 'FontSize', 14);
    % legend("location", "northeastoutside", 'FontSize', 12);   % too many eigenfunctions, they are not labeled
    hold off;
    waitfor(p2)   % Stops the code execution until you close the plot
end


% 5. Displaying values after application of boundary conditions (quantization)
for n = 1:nmax
    kk(n) = n * pi / L;     % Wave number [m^-1] (eq. 2.64)
end
omega = kk * cc;            % Angular frequency [Hz] (eq. 2.65)

% gg study

for n = 1:nmax
    OO(n) = sqrt(omega(n) ^ 2 - gg ^ 2);  % Omega [Hz] (eq. 2.70)
end

fprintf('\nEigenvalues\n');
fprintf('   n        kk(n)[1/m]      OO_n[Hz]        nu_n[Hz]         T_n[s]\n');
for n = 1:nmax
    % Remember that experimental input data is limited to
    % fprintf('%4d %#15.3G %#15.3G %#15.3G %#15.3G\n', n, kk(n), OO(n), OO(n)/(2*pi), 2*pi/OO(n));
    a1(n, 1) = n;
    a1(n, 2) = kk(n);
    a1(n, 3) = OO(n);
    a1(n, 4) = OO(n)/(2*pi);
    a1(n, 5) = 2*pi/OO(n);
end
printtable(a1, 'LaTex', true);

% Omegas study
fprintf('\Omegas study\n');
for gg = 100:300:1000
    gg
    fprintf('   n        oo_n[Hz]      OO_n[Hz]');
    for n = 1:nmax
        b1(n, 1) = n;
        b1(n, 2) = omega(n);
        b1(n, 3) = sqrt(omega(n) ^ 2 - gg ^ 2);
    end
printtable(b1, 'LaTex', true);
end

gg = 100;  % Reset gg to the original value

% 6. TEST ORTHONORMALIZATION OF EIGENFUNCTIONS
fprintf('\nOrthonormalization test\n');
fprintf('   n     <phi_n|phi_n>\n');
for (n=1:nmax)
    % fprintf('%4d %#15.6G\n', n, trapz(x, phi(n, :) .* phi(n, :)));
    a2(n, 1) = n;
    a2(n, 2) = trapz(x, phi(n, :) .* phi(n, :));
end
printtable(a2, 'LaTex', true);

% 7. OVERLAP INTEGRALS BETWEEN INITIAL CONDITIONS AND THE EIGENFUNCTIONS FOR THE VARIOUS PINCH TYPES
pinch = 2  % Pinch type 2 (eq. 2.45), according to the task

fprintf('\nOverlap integrals\n');
fprintf('   n     <phi_n|f>          <phi_n|g>\n');
for (n=1:nmax)
    % fprintf('%4d %#15.6G %#15.6G\n', n, trapz(x, f(pinch, :) .* phi(n, :)), trapz(x, g(pinch, :) .* phi(n, :))); % Overlap integrals (eq. 2.38)
    a3(n, 1) = n;
    a3(n, 2) = trapz(x, f(pinch, :) .* phi(n, :));
    a3(n, 3) = trapz(x, g(pinch, :) .* phi(n, :));
end
printtable(a3, 'LaTex', true);

% 8. COMPARISON OF F AND APPROX F
% Variable declaration
f_approx = 0;
g_approx = 0;

if (strcmp(graph_f_approx,'yes') == 1);
    p3 = plot(x / L, f(pinch, :), ';f(x) [m];', 'LineWidth', 0.5);
    hold on;
    for n=1:nmax
      f_approx = f_approx + trapz(x, f(pinch, :) .* phi(n,:)) * phi(n,:);
    end
    p3 = plot(x / L, f_approx, ';f{_{approx}}(x) [m];', 'LineWidth', 0.5, 'LineStyle', '--');

    p3 = plot(x / L, g(pinch, :), ';g(x) [m/s];', 'LineWidth', 0.5);
    for n=1:nmax
      g_approx = g_approx + trapz(x, g(pinch, :) .* phi(n,:)) * phi(n,:);
    end
    p3 = plot(x / L, g_approx, ';g{_{approx}}(x) [m/s];', 'LineWidth', 0.5, 'LineStyle', '--');

    xlim([-(xmax / L) * 0.1, (xmax / L) * 1.1]);
    ylim([-1.1, 1.1]);
    xlabel('x/L', 'FontSize', 14);
    ylabel('y(x)', 'FontSize', 14);
    title('Comparison of initial conditions and approximations', 'FontSize', 14);
    legend("location", "northeast", 'FontSize', 14);
    hold off;
    waitfor(p3)   % Stops the code execution until you close the plot
end

% 9. Numerical array of the spatial discretization of phi and successive amplitude paterns
% Overlap integrals between initial conditions and the eigenfunctions
for n = 1:nmax
    phi_f(n) = trapz(x, f(pinch, :) .* phi(n,:));
    phi_g(n) = trapz(x, g(pinch, :) .* phi(n,:));
end

% Wave evolution function
if (strcmp(graph_time,'yes') == 1);    
    for i = 1:nst
        for j = 1:npt
            psi_plot(j, i) = 0;
            for n = 1:nmax
                term = phi_f(n) * cos(OO(n) * tt(i)) + (1/OO(n)) * phi_g(n) * sin(OO(n) * tt(i));
                psi_plot(j, i) = psi_plot(j, i) + phi(n, j) * e ^ (-gg * tt(i)) * term;  % Problem solution (eq. 2.82)
            end
        end
    end

    % PLOT OF THE WAVE EVOLUTION
    for (i=1:nst)
        p4 = plot(x/L, psi_plot(:, i), strjoin({';t=', num2str(tt(i)), 's;'}, ' '), 'LineWidth', 0.7);
        hold on;
    end
    xlim([-(xmax / L) * 0.1, (xmax / L) * 1.1]);
    ylim([-1.1, 1.1]);
    xlabel('x/L', 'FontSize', 14);
    ylabel('{\psi}(x, t) [m]', 'FontSize', 14);
    title('Wave evolution frames', 'FontSize', 14);
    legend("location", "northeastoutside", 'FontSize', 14);
    hold off;
    waitfor(p4)     % Stops the code execution until you close the plot
end

% 10. MOVIE
if (strcmp(graph_movie,'yes') == 1);
    nmax = 10       % Maximum number of modes (to be adjusted for the sake of the movie)
    shortest_period = (2 * pi) / max(OO)    % Determine the shortest period
    movie = round((tmax / shortest_period) * 10)    % Number of frames for movie, 10 frames per shortest period
    stepmov = (tmax - tmin)/(movie - 1); % Step t size
    normx = x / L;       % Normalized x discretized abscissa

    for i = 1:movie
        tmov(i) = tmin + (i - 1) * stepmov;
    end

    % Initial psi t = 0 (eq. 2.43)
    for i = 1:npt
        psi(i) = 0;
        for n = 1:nmax
            psi(i) = psi(i) + phi_f(n) * phi(n, i);
        end
    end

    g5 = plot(normx, psi, 'XDataSource', 'normx', 'YDataSource', 'psi', 'linestyle', '-', 'linewidth', 2);

    xlim([-(xmax / L) * 0.1, (xmax / L) * 1.1]);
    ylim([-1.1, 1.1]);
    line([xmin/L,xmax/L], [0, 0], 'linestyle', '-', 'linewidth', 1, 'color', [0.5, 0.5, 0.5]);
    xlabel ('x/L', 'FontSize', 14);
    ylabel ('{\psi}(x, t) [m]', 'FontSize', 14);
    title ('Wave evolution movie', 'FontSize', 14);

    % Loop on frames

    for i = 1:movie
        pause(0);
        for j = 1:npt
            psi(j) = 0;
            for n = 1:nmax
                term = phi_f(n) * cos(OO(n) * tmov(i)) + (1 / OO(n)) * phi_g(n) * sin(OO(n) * tmov(i));
                psi(j) = psi(j) + phi(n, j) * e ^ (-gg * tmov(i)) * term;  % Problem solution (eq. 2.82)
            end
        end
        refreshdata();
    end
end