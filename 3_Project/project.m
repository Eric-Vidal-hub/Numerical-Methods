clear all;  % Clear all variables
clc;        % Clear the command window
close all;  % Close all figures
source('mystartdefaults.m'); % Contains SI physical sonstants

tic

% Units for free electrons in vacuum
recipunit = 1.0E+10;
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm))/qel

% Grid for perturbation
xp_min = 0
xp_max = 1
n=100   % Number of points
step = (xp_max-xp_min)/n

% Grid for wavefunction
x_min = -10
x_max = 10

m=floor((x_max-x_min)/step)

E_0 = 1     % Energy of the particle [eV] (eq. 4.19)
k_0 = sqrt(E_0/ekinscale)   % Wave vector [1/A] (eq. 4.19)
lambda = (2*pi)/k_0        % De Broglie Wavelength [A]

% Solving inside the perturbation
xp = zeros(n,1);    % Fixiing the size of the array as a column vector
for i=1:n        % Remember: xp_min and x_max must not be included
    xp(i) = xp_min + (i-1)*step + step/2; % Discretized abscissa [A]
end

Phi0p = exp(1i * k_0 * xp); % Initial wavefunction inside the perturbation (eq. 4.23)

U=zeros(n,1);
for i=1:n
    U(i) = 2; % Discretization inside the perturbation [eV] (eq. 4.9)
end

V = U/ekinscale; % Unit conversion in k^2 [1/A^2] (eq. 4.9)

G0 = zeros(n,n); % Green's function matrix
for i=1:n
    for j=1:n
        G0(i,j) = step * exp(1i * k_0 * abs(xp(i)-xp(j))) / (2 * 1i * k_0);
    end
end

T = eye(n,n) - G0 * diag(V); % Matrix in eq. 4.51

Phip = T \ Phi0p; % Wavefunction inside the perturbation (eq. 4.51)

for i=1:n
    Phip2(i) = abs(Phip(i)^2); % Wavefunction inside the perturbation (eq. 4.51)
end

% Check graph of result inside perturbation
hf(1) = figure('NumberTitle', 'off', 'name', ['Square modulus of the wavefunction inside the perturbation'])
plot(xp,Phip2); % Plotting the wavefunction inside the perturbation
xlabel('x [A]'); % Label for the x axis
ylabel('|\Phi|^2 [1/A]'); % Label for the y axis    
movegui(hf(1),'west') % Move the figure to the northwest of the screen

% Propagating the solution outside the perturbation
% does not require large matrices (eq. 4.49 and 4.36)
x = zeros(m,1); Phis = zeros(m,1); Phi = zeros(m,1); Probab = zeros(m,1);

for i = 1:m
    x(i) = x_min + (i-1)*step + step/2; % Discretized abscissa [A], x = xp inside perturbation
    Phis(i) = 0; % Initial wavefunction outside the perturbation
    for j=1:n
        Phis(i) = Phis(i) + step * exp(1i * k_0 * abs((x(i) - xp(j))))/(2 * 1i * k_0) * V(j) * Phip(j); % Wavefunction outside the perturbation (eq. 4.49)
    end
    Phi(i) = exp(1i * k_0 * x(i)) + Phis(i); % Total wavefunction (eq. 4.36)
    Proba(i) = abs(Phi(i)^2); % Probability density
end

Ref = abs(Phis(1)/exp(1i * k_0 * x(1)))^2; % Reflection coefficient
Tra = Proba(m)
RpT = Ref/Tra % Ratio between reflection and transmission coefficients
toc

% For plotting: Curve of potential in [x_min, x_max]
Ugraph = zeros(m,1);
for i=1:m
    if (x(i) >= xp_min) && (x(i) <= xp_max)
        for j=1:n
            if (abs(x(i) - xp(j)) < 1e-12)
                Ugraph(i) = U(j); % Potential in [x_min, x_max] [eV]
            end
        end
    end
end
    % Check graph of result outside perturbation
hf(2) = figure('NumberTitle', 'off', 'name', ['Scattering eigenstate']);
hold on
plot(x,Proba); % Plotting the wavefunction outside the perturbation
plot(x,Ugraph); % Plotting the potential in [x_min, x_max] [eV]
plot(x, real(Phi)); % Plotting the real part of the wavefunction outside the perturbation
plot(x, imag(Phi)); % Plotting the imaginary part of the wavefunction outside the perturbation
xlabel('x [A]'); % Label for the x axis
legendstr{1} = ['|\Phi(x)|^2']; % Legend for the wavefunction
legendstr{2} = ['U(x) [eV]']; % Legend for the potential
legendstr{3} = ['Re(\Phi(x))']; % Legend for the real part of the wavefunction 
legendstr{4} = ['Im(\Phi(x))']; % Legend for the imaginary part of the wavefunction
lgd = legend(legendstr, 'location', 'eastoutside'); % Legend for the wavefunction
movegui(hf(2),'east') % Move the figure to the northeast of the screen


% Current of probability in QM
% J(r,t) = (hbar/(2*m*1i)) * (conj(Phi) * grad(Phi) - Phi * grad(conj(Phi)))
% J(r,t) = (hbar/m) * real(conj(Phi) * grad(Phi)/1i)

Current = zeros(n,1);
deriv = zeros(n,1);
for i=1:m-1
    deriv(i) = (Phi(i+1) - Phi(i))/step; % Derivative of the wavefunction on the discretized grid
end

deriv(m) = deriv(m-1); % Derivative of the wavefunction on the discretized grid

speed_factor = recipunit * hbar/elm; % Conversion factor for the speed of the particle
speed_factor = speed_factor * 1e-6; % Conversion factor for the speed of the particle