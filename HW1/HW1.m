% Task: NM-HW1 2.9.5 Damped Vibrating String
% Deadline: October 4th, 2023 at 12:00.
% Author: Eric Vidal Marcos
% Group: 2 (QuanTEEM)
% Email: ericvidalmarcos@gmail.com
% Supervisors: Dr. Alain Dereux & Dr. Matthieu Sala

% Clear workspace
clear all   % Clear all variables
close all   % Close all variables
clc         % Clear command window

% Plot commands
graph_pinch = 'no';             % Plot initial conditions of f and g
graph_eigenfunctions = 'no';    % Plot the eigenfunctions
% graph_f = 'no';                 % Plot the comparison of f and approx f
% graph_time = 'no';              % Plot the evolution of the wave
% graph_movie = 'no';             % Plot the movie of the wave evolution

% Pinch commands
minpinch = 1;                   % Minimum number of pinches
maxpinch = 2;                   % Maximum number of pinches
% If you want to work with just one pinch, set minpinch = maxpinch

% 1. Input data and discretization parameters
L = 0.328           % String length [m]
mu = 0.660E-3       % String mass per unit length (mu = M / L) [kg/m]
F = 55.0            % Tension force [N]
cc = sqrt(F/mu)     % Speed of wave propagation [m/s] (eq. 2.48)
nmax = 20           % Maximum number of modes (to be adjusted for the good convergence of the approximation)
npt = 201           % Number of sampling points used in the approximation
nst = 13            % Number of sample times
movie=200           % Number of frames for movie

gg = 100.0                  % Gamma coefficient [Hz]
dampf = 2 * gg / (cc ^ 2)   % Damping factor [s/m] (eq. 2.86)
tmax = (2 / gg) * log(10)   % Maximum time [s] (my pdf eq.)
tmin = 0;                   % Minimum time [s]
tstep = (tmax - tmin) / (nst - 1);  % Time step [s]
for i=1:nst
  t = tmin + (i - 1) * tstep;    % Discretized time [s]
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
    elseif (pinch == 2);  % eq. (2.44)
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
        xlabel('x/L');
        ylabel('y(x)');
        title(strjoin({'Pinch', num2str(pinch)}, ' '));
        legend("location", "northeast");
        waitfor(p1)   % Stops the code execution until you close the plot
    end
endif

% 4. Compute the eigenfunctions and plot
for n = 1:nmax
    for i = 1:npt
        phi(n, i) = sqrt(2 / L) * sin(n * pi * x(i) / L);   % Eigenfunctions (eq. 2.75)
    end
end

% GRAPH OF EIGENFUNCTIONS
if (strcmp(graph_eigenfunctions, 'yes') == 1);
    for n = 1:nmax
        p2 = plot(x/L, phi(n, :), 'LineWidth', 1);
        hold on;
    end
    xlim([xmin/L,xmax/L]);
    ylim([-sqrt(2/L),sqrt(2/L)]);
    xlabel('x/L');
    ylabel('{\Phi(x)}');
    title('Eigenfunctions');
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
fprintf('   n        kk(n)[1/m]   OO_n[Hz]     nu_n[Hz]     T_n[s]\n');
for n = 1:nmax
    % Remember that experimental input data is limited to
    fprintf('%4d %#15.3G %#15.3G %#15.3G %#15.3G\n', n, kk(n), OO(n), OO(n)/(2*pi), 2*pi/OO(n));
end