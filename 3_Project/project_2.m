clear all;  % Clear all variables
clc;        % Clear the command window
close all;  % Close all figures
source('mystartdefaults.m');    % Contains SI physical sonstants

tic

% Units for free electrons in vacuum
recipunit = 1.0E+10;
ekinscale = ((hbar * recipunit)^2 / (2.0 * elm))/qel

% Grid for perturbation
xp_min = 0
xp_max = 80
n=100   % Number of points
step = (xp_max-xp_min)/n

E_min = 0
E_max = 10
E_step = 0.1

E_0 = [E_min:E_step:E_max];
k_0 = sqrt(E_0 / ekinscale);

% Solving inside the perturbation
xp = zeros(n,1);    % Fixiing the size of the array as a column vector
for i=1:n           % Remember: xp_min and x_max must not be included
    xp(i) = xp_min + (i-1)*step + step/2;   % Discretized abscissa [A]
end

U=zeros(n,1);
for i=1:n
    if (xp(i)>0 && xp(i)<15)
        U(i) = 0.1;
    else
        U(i) = 0;
    end
end

V = U/ekinscale; % Unit conversion in k^2 [1/A^2] (eq. 4.9)

%% SOLVING INSIDE PERTURBATION

m=2; x=zeros(m,1); Phis=zeros(m,1); Phi=zeros(m,1);

x(1) = xp_min - 1;
x(2) = xp_max + 1;

%% LOOP OVER INCIDENT ENERGIES
G0 = zeros(n,n);  % Greens function matrix inside the barrier

fprintf('    k       E_0         R        T        R+T\n')
for k=1:length(E_0)
    if (E_0(k) == 0)
        Ref(k) = 1;
        Tra(k) = 0;
        RpT(k) = 1;
    else
        Phi0p = exp(1i * k_0(k) * xp);   % Incident plane wave inside the perturbation (eq. 4.23)
        
        for j=1:n
            for i=1:n
                G0(i,j) = step * exp(1i * k_0(k) * abs(xp(i) - xp(j))) / (2 * 1i * k_0(k)); %  (eq. 4.34)
            end
        end

        T = eye(n,n) - G0 * diag(V); % Matrix in eq. (4.51)
        
        Phip=T\Phi0p;
        for i=1:m
            Phis(i) = 0;
            for j=1:n
                Phis(i) = Phis(i) + step * (exp(1i * k_0(k) * abs(x(i) - xp(j))) / (2 * 1i * k_0(k))) * V(j) * Phip(j);
            end
            Phi(i) = exp(1i * k_0(k) * x(i) + Phis(i));
        end
        Ref(k) = abs(Phis(1)/exp(1i * k_0(k) * x(1)))^2;     % Reflection coefficient
        Tra(k) = abs(Phi(2))^2;                             % Transmission coefficient    
        RpT(k) = Ref(k) + Tra(k);                           % Sanity check R+T=1
    end
    fprintf('%5d  %15.6G  %15.6G  %15.6G  %15.6G\n',k,E_0(k),Ref(k),Tra(k),RpT(k))
end

toc

%% Plot
