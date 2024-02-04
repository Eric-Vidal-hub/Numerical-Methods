% Propogate solution out from peturbation to entire grid
%% Solution outside peturbation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phis = zeros(L1,1);
phitot = zeros(L1,1);

for i = 1:L1 % Equation 4.13 - Lippmann-Schwinger equation, identical to Huygen's principle
    phis(i) = 0;
    for j = 1:L2 % integrating Green's function
        phis(i) = phis(i) + (xres * (Greenfunc(x(i),Ux(j),k0,k1) * (W(j)) * phip(j)));
    end
    if(x(i)>=0) % adding incicdent/transmitted field
        phitot(i) = tb*exp(1i*k1*x(i)) + phis(i);                       % total on RHS = transmitted from step + scattered
    else
        phitot(i) = (rb)*exp(-1i*k0*x(i)) + exp(1i*k0*x(i)) + phis(i);
    end
end

if(plotfullwavefn=='y')
    fig2 = figure;
    hold on
    grid on
    plot(x,abs(phitot).^2)
    plot(x,(Vfull))
    titlestr = sprintf('Solution of scattering states from %d to %d (Ang)',xmin, xmax);
    title(titlestr,'FontSize',20)
    xlabel('x (angstrom)','FontSize',15);
    ylabel('abs(psi)^2','FontSize',15);
end

%% transmission and reflection coefficients
Trans = (real(k1)/real(k0)) * abs(phitot(end))^2;
Refl = abs((rb*exp(-1i*k0*x(1)) + phis(1)) / exp(1i*k0*x(1)))^2;

fprintf('\n Reflection coefficient = %d \n Transmission coefficient = %d\n Sum of reflection and transmission coefficients is %15.5e\n',Refl, Trans,(Refl+Trans))
