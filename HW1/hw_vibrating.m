clear all;  % Clear all variables
clc;        % Clear the command window

% Plot commands
graph_pinch = 'no';             % Plot initial conditions of f and g
graph_eigenfunctions = 'no';    % Plot the eigenfunctions
graph_f = 'no';                 % Plot the comparison of f and approx f
graph_time = 'no';              % Plot the evolution of the wave
graph_movie = 'yes';            % Plot the movie of the wave evolution

% 1. Define the input data
L = 0.328           % String length [m]
mu = 0.660E-3       % String mass per unit length [kg/m]
F = 55.0            % Tension force [N]
cc = sqrt(F/mu)     % Speed of wave propagation [m/s] (eq. 2.2)
pinch = 1           % Pinched point (0 = no pinch, 1 = pinch at the center)
nmax = 20           % Number of modes to be considered, eigenfunctions used in the approximation
npt = 201           % Number of sampling points used in the approximation
nst = 13            % Number of time samples for the longest period of the spectrum
movie = 200         % Number of frames in the movie

% 2. Array of discretized value of x
xmin = 0;       % Minimum abscissa
xmax = L;       % Maximum abscissa
step = (xmax - xmin)/(npt - 1); % Step size
for i = 1:npt
    x(i) = xmin + (i - 1) * step; % Descretized abscissa
end 

% 3. Compute numerical arrays for both pinches and to be used as initial conditions
if (pinch == 0)
    fprintf('Pinch point not defined\n');
    return;
end

if (pinch==1);
  fprintf('The pinch 1 has been chosen.')
    p=L/2;   % Pinched point
    z = 1;   % Amplitude of the wave
    for (i=1:npt)
        if (x(i) <= p)
            f(i) = x(i) * z / p;
        else
            f(i) = z * (x(i) - L) / (p-L);
        end
    end
    for i=1:npt     % Initial distribution of velocity field
        g(i) = 0;
    end
end

if (pinch==2);
  fprintf('The pinch 2 has been chosen.');
    p1 = L/4.;    % Pinched point 1
    z1 = 1;       % Amplitude of the wave 1

    p2 = (3./4.) * L; % Pinched point 2
    z2 = -z1;       % Amplitude of the wave 2
    for (i=1:npt)
        if (xmin <= x(i)) && (x(i) <= p1)
            f(i) = x(i) * z1 / p1;
        elseif (p1 < x(i)) && (x(i) <= p2)
            f(i) = (z1 - z2) / (p1 - p2) * (x(i) - p1) +z1;
        elseif (p2 < x(i)) && (x(i) <= xmax)
            f(i) = z2 * (x(i) - L) / (p2 - L);
        end
    end
    for i=1:npt     % Initial distribution of velocity field
        g(i) = 0;
    end
end

if (pinch==3);
    m = 3   % Proportional constant
    for i=1:npt % Deformation = an aigenmode
        f(i) = 1.5 * sqrt(2/L) * sin(m*pi*x(i)/L);
    end
end

if (strcmp(graph_pinch,'yes') == 1);
  p1 = plot(x/L, f);
  hold on
  p1 = plot(x/L, g);
  hold off
  xlabel('x/L');
  ylabel('f(x) and g(x)');
  title('Pinch');
  waitfor(p1)   % Stops the code execution until you close the plot
endif

% 4. EIGENFUNCTIONS
for n=1:nmax
    for i=1:npt
        phi(n, i) = sqrt(2/L) * sin(n*pi*x(i)/L);
    end
end

% GRAPH OF EIGENFUNCTIONS
if (strcmp(graph_eigenfunctions,'yes') == 1);
    for (n=1:nmax)
        p2 = plot(x/L, phi(n,:), 'LineWidth', 2);
        hold on;
    end
    xlim([xmin/L,xmax/L]);
    ylim([-sqrt(2/L),sqrt(2/L)]);
    xlabel('x/L');
    ylabel('phi(x)');
    title('Eigenfunctions');
    hold off;
    waitfor(p2)   % Stops the code execution until you close the plot
end

% 5. EIGENVALUES
for n=1:nmax
    vk(n)=n*pi/L;
end
omega = vk * cc;
fprintf('\nEigenvalues\n');
fprintf('   n        vk(n)[1/m]   omega(n)[Hz]     nu(n)[Hz]     period(n)[s]\n');
for n=1:nmax
    % Remember that experimental input data is limited to
    fprintf('%4d %#15.3G %#15.3G %#15.3G %#15.3G\n', n, vk(n), omega(n), omega(n)/(2*pi), 2*pi/omega(n));
end

% 6. TEST ORTHONORMALIZATION OF EIGENFUNCTIONS

fprintf('\nOrthonormalization test\n');
fprintf('   n   <phi_n|phi_n>\n');
for (n=1:nmax)
    fprintf('%4d %#15.6G\n', n, trapz(x, phi(n,:) .* phi(n,:)))
end

% 7. OVERLAP INTEGRALS BETWEEN INITIAL CONDITIONS AND THE EIGENFUNCTIONS FOR THE VARIOUS PINCH TYPES
fprintf('\nOverlap integrals\n');
fprintf('   n   <phi_n|f>    <phi_n|g>\n');
for (n=1:nmax)
    fprintf('%4d %#15.6G %#15.6G\n', n, trapz(x, f .* phi(n,:)), trapz(x, g .* phi(n,:)));
end

% 8. COMPARISON OF F AND APPROX F

% Variable declaration
f_approx = 0;
g_approx = 0;

if (strcmp(graph_f,'yes') == 1);
    p3 = plot(x, f, 'LineWidth', 0.5);
    hold on;
    for n=1:nmax
      f_approx = f_approx + trapz(x, f .* phi(n,:)) * phi(n,:);
    endfor
    p3 = plot(x,f_approx, 'LineWidth', 0.5, 'LineStyle', '--');

    p3 = plot(x, g, 'LineWidth', 0.5);
    for n=1:nmax
      g_approx = g_approx + trapz(x, g .* phi(n,:)) * phi(n,:);
    endfor
    p3 = plot(x,g_approx, 'LineWidth', 0.5, 'LineStyle', '--');

    xlabel('x');
    ylabel('fs');
    title('Comparison of fs and gs');
    hold off;
    waitfor(p3)   % Stops the code execution until you close the plot
end

% 9. Numerical array of the spatial discretization of phi and successive amplitude paterns
longest_period = (2 * pi) / min(omega); % Longest period of the spectrum
tmin = 0;                   % Minimum t abscissa
tmax = longest_period;      % Maximum t abscissa
step = (tmax - tmin)/(nst - 1); % Step t size

% Overlap integrals between initial conditions and the eigenfunctions for the various pinch types
for n = 1:nmax
    phi_f(n) = trapz(x, f .* phi(n,:));
    phi_g(n) = trapz(x, g .* phi(n,:));
end

if (strcmp(graph_time,'yes') == 1);    
    for i = 1:nst
        t = tmin + (i - 1) * step; % Descretized t abscissa
        for j = 1:npt
            psi(j, i) = 0;
            for n = 1:nmax
                buf = phi_f(n) * cos(omega(n) * t) + (1/omega(n)) * phi_g(n) * sin(omega(n) * t);
                psi(j, i) = psi(j, i) .+ phi(n, j) * buf;
            endfor
        endfor
    end

    % PLOT STARTS

    for (i=1:nst)
        p4 = plot(x/L, psi(:, i), 'LineWidth', 2);
        hold on;
    end
    % xlim([xmin/L,xmax/L]);
    % ylim([-sqrt(2/L),sqrt(2/L)]);
    xlabel('x/L');
    ylabel('psi(x, t)');
    title('Evolution');
    hold off;
    waitfor(p4)   % Stops the code execution until you close the plot
end

% 10. MOVIE
if (strcmp(graph_movie,'yes') == 1);
    longest_period = (2 * pi) / min(omega); % Longest period of the spectrum
    tmin = 0;                       % Minimum t abscissa
    tmax = 2 *longest_period;       % Maximum t abscissa
    step = (tmax - tmin)/(movie - 1); % Step t size
    normx = x / L;       % Normalized x discretized abscissa

    % Initial psi t = 0 [eq. 2.43]
    for i = 1:npt
        psi(i) = 0;
        for n = 1:nmax
            psi(i) = psi(i) + phi_f(n) * phi(n, i);
        end
    end

    g5 = plot(normx, psi, 'XDataSource', 'normx', 'YDataSource', 'psi', 'linestyle', '-', 'linewidth', 2);

    xlim([xmin/L,xmax/L]);
    ylim([-sqrt(2/L),sqrt(2/L)]);
    line([xmin/L,xmax/L], [0, 0], 'linestyle', '-', 'linewidth', 1, 'color', [0.5, 0.5, 0.5]);
    xlabel ('x/L');
    ylabel ('Amplitude');
    title ('Wave evolution');

    % Loop on frames

    for k = 1:movie
        pause(0);           % Wait for a key press
        t = tmin + (k - 1) * step;  % Descretized t abscissa
        for j = 1:npt
            psi(j) = 0;
            for n = 1:nmax
                buf = phi_f(n) * cos(omega(n) * t) + (1/omega(n)) * phi_g(n) * sin(omega(n) * t);
                psi(j) = psi(j) .+ phi(n, j) * buf;
            endfor
        end
        refreshdata();
    end;
end