%% Calculation of Reflection-Absorption-Transmission 
% single script allows to vary incident energy as well as 
% peturbation and bias
% still pending some sanity checks on tuning reference
clearvars
tic
%% Plot flags
plotpotential = 'n  ';
plotpeturbsoln = 'n';
plotfullwavefn = 'n';

%% select RAT curve
makespectrum = 'n';
tunereference = 'y';
tunestep1 = 'n';
tunestep2 = 'n';

% range of whicever parameter is being varied
Estart = -0.2;
Eend = 0.2;
Eres = 0.001;  % resolution of spectrum

%% default values
% Will be used for parameters not being tuned
mystartdefaults
recipunit = 1e10;          % Reciprocal space unit = 1e10 1/m
% Constant to convert k (A^2) to Energy (eV)
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;
fprintf('\n ekinscale = %d \n',ekinscale);
spatial_disc

defaultE0 = 0.01;

UL = 0;       % eV Potential of reference system if x<0
UR = -0.1;    % eV Potential of reference system if ><0

U1 = 0.2;       % eV                   % Height of 1st barrier
U2 = 0.2;       % eV                   % Height of 2nd barrier\
peturb_init
W = W/ekinscale;

lifetime = 1e-9;
gamma = hbar*(2*pi/lifetime)/qel;
fprintf('\n Gamma  = %d\n',gamma)
%% switch case to establish variables for tuning
if(makespectrum=='y' && not(tunereference=='y') && not(tunestep1=='y') && not(tunestep2=='y'))
    E0range = (Estart+Eres/2:Eres:Eend-Eres/2);
    L3 = length(E0range);
    fprintf('\n Calculating Energy Spectrum from Estart = %d eV to Eend = %2.2d eV; step = %d \n',Estart,Eend, UR-UL)
    plt_title = sprintf('Energy spectrum in range [%2d,%2d]; step = %d',Estart,Eend, UR-UL);
    xvar = E0range;
elseif(not(makespectrum=='y') && (tunereference=='y') && not(tunestep1=='y') && not(tunestep2=='y'))
    URrange = (Estart+Eres/2:Eres:Eend-Eres/2);
    L3 = length(URrange);
    fprintf('\nCalculating R/A/T for tuning Bias from Estart = %d eV to Eend = %d eV\n',Estart,Eend)
    plotpotential ='n';     % necessary
    plt_title = sprintf('Tuning bias in range [%2d,%2d]; E0 = %d',Estart,Eend,defaultE0);
    xvar = -URrange;
elseif(not(makespectrum=='y') && not(tunereference=='y') && (tunestep1=='y') && not(tunestep2=='y'))
    U1range = (Estart+Eres/2:Eres:Eend-Eres/2);
    L3 = length(U1range);
    fprintf('\nCalculating R/A/T for tuning step 1\n')
    plotpotential ='n';     % necessary
    plt_title = sprintf('Tuning first step in range [%2d,%2d]; E0 = %d',Estart,Eend, defaultE0);
    xvar = U1range;
elseif(not(makespectrum=='y') && not(tunereference=='y') && not(tunestep1=='y') && (tunestep2=='y'))
    U2range = (Estart+Eres/2:Eres:Eend-Eres/2);
    L3 = length(U2range);
    fprintf('\nCalculating R/A/T for tuning step 2 from Estart = %d eV to Eend = %d eV\n',Estart,Eend)
    plotpotential ='n';     % necessary
    plt_title = sprintf('Tuning second step in range [%2d,%2d]; E0 = %d',Estart,Eend, defaultE0);
    xvar = U2range;
else
    warning('Terminating Program: Please select exactly one parameter for R/A/T curve')
    return
end

%% Preallocation
% significantly boosts performance
trans=zeros(L3,1);
refl=zeros(L3,1);
absorb=zeros(L3,1);

%% Calculation of Spectrum

for n = 1:L3
    if(makespectrum=='y')
        E0 = E0range(n);
    end
    if(tunereference=='y')
        E0 = defaultE0;
        UR = URrange(n);
        peturb_init
        W = W/ekinscale;
    end
    if(tunestep1=='y')
        E0 = defaultE0;
        U1 = U1range(n);
        peturb_init
        W = W/ekinscale;
    end
    if(tunestep2=='y')
        E0 = defaultE0;
        U2 = U2range(n);
        peturb_init
        W = W/ekinscale;
    end
    k0 = sqrt((E0-UL+(1i*gamma))/ ekinscale); % Eq. A.58
    k1 = sqrt((E0-UR+(1i*gamma))/ ekinscale); % Eq. A.58

    rb = (k0-k1)/(k1+k0); % Eq. A.64
    tb = (2*k0)/(k1+k0); % Eq. A.64

    % Incident wavefn
    phi0=tb*exp(1i*k1*Ux); % Eq. A.65

    %Green's function
    G = Greenfunc(Ux,Ux,k0,k1); 

    % Solve System - Eq 4.51
    T = eye(L2,L2) - xres*G*diag(W); 
    phip = T\phi0;   
    % calc_peturbsoln

    % Huygens principle
    phis = zeros(2,1);
    phitot = zeros(2,1);
    testp = [x(1),x(end)];      % working with just 2 points because that is all you need for the spectrum
    
    for pp = 1:2 % Equation 4.13 - Lippmann-Schwinger equation, identical to Huygen's principle
        phis(pp) = 0;
        for qq = 1:L2
            phis(pp) = phis(pp) + xres*Greenfunc(testp(pp),Ux(qq),k0,k1) * (W(qq)) * phip(qq);
        end
        if(testp(pp)>0)
            phitot(pp) = tb*exp(1i*k1*testp(pp)) + phis(pp); 
            % total on RHS = transmitted from step + scattered
        else
            phitot(pp) = (rb)*exp(-1i*k0*testp(pp)) + exp(1i*k0*testp(pp)) + phis(pp);  
            % total on LHS = reflected from step + incident + scattered
        end
    end

    trans(n) = (real(k1)/real(k0)) * abs(phitot(2))^2;
    refl(n) = abs((rb*exp(-1i*k0*testp(1)) + phis(1)) / exp(1i*k0*testp(1)))^2;
    absorb(n) = 1 - (refl(n)+trans(n));
end

%% Plots
figure
hold on
plot(xvar,abs(trans),'LineWidth',1.5)
plot(xvar,abs(refl),'LineWidth',1.5)
plot(xvar,abs(absorb),'LineWidth',1.5)
grid on
xlabel('Energy (eV)','FontSize',15)
legend('trans','refl','abs','FontSize',10)
title(plt_title,'FontSize',15);
toc