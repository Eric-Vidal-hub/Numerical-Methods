%% For solving peturbation and plotting wavefunction across sampling grid
% solves wavefunction within peturbation as well as propogating to rest of
% spatial grid
% Calls spatial_grid, peturb_init and mystartdefaults
clearvars
%% Plot flags
plotpotential = 'y';
plotpeturbsoln = 'y';
plotfullwavefn = 'y';

%% units and constants
mystartdefaults
recipunit = 1e10;          % Reciprocal space unit = 1e10 1/m
% Constant to convert k (A^2) to Energy (eV)
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;
fprintf('\n ekinscale = %d \n',ekinscale);

%% Discretisation
spatial_disc

%% define peturbation
UL = 0;       % eV Potential of reference system if x<0
UR = -0.1;    % eV Potential of reference system if ><0

U1 = 0.2;       % eV                   % Height of 1st barrier
U2 = 0.2;       % eV                   % Height of 2nd barrier
peturb_init
%% solve wavefunction within peturbation
E0 = 0.1;             % eV
lifetime = 1e-9;
gamma = hbar*(2*pi/lifetime)/qel;

calc_peturbsoln

%% propogation of wavefunction to complete space
calc_propogation