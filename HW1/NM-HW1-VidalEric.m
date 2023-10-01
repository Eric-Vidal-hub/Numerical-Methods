% Task: 2.9.5 HW1 Damped Vibrating String
% Deadline: October 4th, 2023 at 12:00.
% Author: Eric Vidal Marcos

% Clear workspace
clear all
close all
clc

% 1. Input data and discretization parameters
L = 0.328           % String length [m]
mu = 0.660E-3       % String mass per unit length (mu = M / L) [kg/m]
F = 55.0            % Tension force [N]
cc = sqrt(F/mu)     % Speed of wave propagation [m/s] (eq. 2.48)
gg = 100.0          % Gamma coefficient [Hz]
dampf = 2 * gg / (cc ^ 2)   % Damping factor [s/m] (eq. 2.86)
tmax = (2 / gg) * log(10)   % Maximum time [s] (my pdf eq.)
tmin = 0;
tstep = (tmax - tmin) / (nst - 1);

for i=1:nst
  t = tmin + (k  -1) * timestep;
end

nmax = 20           % Maximum number of modes (to be adjusted for the good convergence of the approximation)
npt = 201           % Number of sampling points used in the approximation
nst = 13            % Number of sample times
movie=200           % Number of frames for movie

% 2. Array of discretized value of x
xmin = 0;           % Minimum abscissa  
xmax = L;           % Maximum abscissa  
xstep = (xmax - xmin) / (npt - 1);
for i=1:npt
  x(i) = xmin + (i-1) * xstep ; % Discretized abscissas  
end

% 3. Compute numerical arrays for both pinches and to be used as initial conditions

for pinch=1:2
    if (pinch==1);
        p=L/2;   % Pinched point
        z = 1;   % Amplitude of the wave
        for (i=1:npt)
            if (x(i) <= p)
                f(pinch, i) = x(i) * z / p;
            else
                f(pinch, i) = z * (x(i) - L) / (p-L);
            end
        end
        for i=1:npt     % Initial distribution of velocity field
            g(pinch, i) = 0;
        end
    elseif (pinch==2);
        p1 = L/4.;    % Pinched point 1
        z1 = 1;       % Amplitude of the wave 1

        p2 = (3./4.) * L; % Pinched point 2
        z2 = -z1;       % Amplitude of the wave 2
        for (i=1:npt)
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