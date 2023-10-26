clear all;  % Clear all variables
clc;        % Clear the command window
close all;  % Close all figures
source('mystartdefaults.m'); % Contains SI physical sonstants

tic

% Units for free electrons in vacuum
recipunit = 1.0E+10;
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm))/qel;

% Grid for perturbation
xp_min = 0;
xp_max = 1;
n=100;   % Number of points
step = (xp_max-xp_min)/n;

% Grid for wavefunction
x_min = -10;
x_max = 10;

m=floor((x_max-x_min)/step);

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
    U(i) = 2; % Discretixation inside the perturbation [eV] (eq. 4.9)
end

V = U/ekinscale; % Unit conversion in k^2 [1/A^2] (eq. 4.9)

G0 = zeros(n,n); % Green's function matrix
for i=1:n
    for j=1:n
        G0 = step * exp(1i * k_0 * abs(xp(i)-xp(j))) / (2 * 1i * k_0);
    end
end

T = eye(n,n) - G0 * diag(V); % Matrix in eq. 4.51

Phip = T \ Phi0p; % Wavefunction inside the perturbation (eq. 4.51)

plot(xp,abs(Phip).^2,'r'); % Plotting the wavefunction inside the perturbation

toc