#! /opt/local/bin/octave
% source('~/bin/mystartdefaults.m');

% vibstr.m : 1-DIM IDEAL (UNDAMPED & PERFECTLY ELASTIC) VIBRATING STRING

			       % OVERWRITE PLOT CONFIGURATION DEFAULTS

set(groot,'defaultAxesXlabel','x/L');
set(groot,'defaultAxesYlabel','Amplitude');

				% GRAPHICAL OUTPUT COMMANDS
graph_eigenfunctions='yes';
graph_pinch='yes';
graph_approx='yes';
graph_sample_times='yes';
graph_movie='yes';

				% INPUT DATA 

L =  0.328       % String length [m]
mu = 0.660E-3    % String mass per unit length [kg/m]
mass = mu*L      % String mass [kg]
F =   55.0       % Tension force applied on string [N]  
c =  sqrt(F/mu)  % Speed [m/s] eq. (2.2)
pinch = 1        % Type of pinch
nmax=20          % Number of eigenfunctions  used in the approximation
npt = 201        % Number of sampling points used in the approximation
nst = 7          % Number of sample times
movie=200        % Number of frames for movie

xmin = 0;        % Minimum abscissa  
xmax = L;        % Maximum abscissa  
step = (xmax-xmin)/(npt-1);
for i=1:npt
  x(i) = xmin + (i-1) * step ; % Discretized abscissas  
end

				% INITIAL CONDITIONS
pinch_defined=0;

if (pinch==1); % eq. (2.44)
  z=1;
  p=L/2;
  for i=1:npt % Deformation (Pinch) at t=0
    if (x(i)<=p);
      f(i) = z*x(i)/p;
    else
      f(i) = z*(x(i)-L)/(p-L);
    end
  end
  
  for i=1:npt % Initial distribution of velocity field
    g(i)=0;
  end
  pinch_defined=1;
end

if (pinch==2); % eq. (2.45)
  z1=1;
  z2=-z1;
  p1=L/4;
  p2=3*L/4;
  for i=1:npt % Deformation (Pinch) at t=0
    if (x(i) <= p1);
      f(i) = z1*x(i)/p1;
    elseif ( p1 < x(i) && x(i) <= p2);
      f(i) =((z1-z2)/(p1-p2)) *(x(i)-p1) + z1 ;
    elseif ( x(i) > p2);
      f(i) = z2*(x(i)-L)/(p2-L);
    end
  end
  
  for i=1:npt % Initial distribution of velocity field
    g(i)=0;
  end
  z=max(abs(z1),abs(z2));
  pinch_defined=1;
end

if (pinch==3);
  m=3;
  for i=1:npt % Deformation = an eigenmode
    f(i)=sqrt(2/L) * sin(m*pi*x(i)/L); 
  end
  
  for i=1:npt % Initial distribution of velocity field
    g(i)=0;
  end
  z=sqrt(2/L);
  pinch_defined=1;
end

if (pinch_defined==0)
  fprintf('Stop because pinch is not defined\n');
  return;
end
				% EIGENFUNCTIONS

for n=1:nmax
  for i=1:npt
    phi(n,i)= sqrt(2/L) * sin(n*pi*x(i)/L); 
  end
end

				% EIGENVALUES
for n=1:nmax
  vk(n)=n*pi/L;
end
omega=vk*c;
fprintf('\nEigenvalues\n');
fprintf('   n      vk(n)[1/m]     omega(n)[Hz]      nu(n)[Hz]    period(n)[s]\n');
for n=1:nmax
   % Remember that experimental input data are limited to 3 significant digits
   fprintf('%4d %#15.3G %#15.3G %#15.3G %#15.3G\n',n,vk(n),omega(n),omega(n)/(2*pi),2*pi/omega(n));
end
Longest_period  = (2*pi)/(min(omega)); % Longest  period contained in Fourier spectrum
Shortest_period = (2*pi)/(max(omega)); % Shortest period contained in Fourier spectrum

				% ORTHONORMALISATION
for n=1:nmax
  for i = 1:npt
    term(i)= phi(n,i)*phi(n,i);
  end
  phi2(n)=trapz(x,term);     % Trapezoidal integration rule
end
fprintf('\nOrthonormalisation\n');
fprintf('   n      <phi(n)|phi(n)>\n');
for n=1:nmax
  fprintf('%4d %#15.6G \n',n,phi2(n));
end
				% OVERLAP INTEGRALS
for n=1:nmax
  for i = 1:npt
    term(i)= phi(n,i)*f(i);
  end
  phidotf(n)=trapz(x,term);      
  for i = 1:npt
    term(i)= phi(n,i)*g(i);
  end
  phidotg(n)=trapz(x,term); 
end
fprintf('\nOverlap integrals\n');
fprintf('   n      <phi(n)|f>      <phi(n)|g>\n');
for n=1:nmax
  fprintf('%4d %#15.6G %#15.6G\n',n,phidotf(n),phidotg(n));
end
fprintf('\n')

				% HOW GOOD IS THE APPROXIMATION ?
for i=1:npt
  fapprox(i)=0;
  gapprox(i)=0;
  for n=1:nmax
    fapprox(i)=fapprox(i)+phidotf(n)*phi(n,i);
    gapprox(i)=gapprox(i)+phidotg(n)*phi(n,i);
  end
end


				% SOLUTION AT SAMPLE TIMES WITHIN A PERIOD

Longest_period = (2*pi)/(min(omega)); % Longest period contained in Fourier spectrum
tmin = 0;
tmax = Longest_period/2;
timestep = (tmax-tmin)/(nst-1);

for k=1:nst
  t=tmin+(k-1)*timestep;
  for i=1:npt
    Psi(k,i)=0;
    for n=1:nmax
      buf = phidotf(n)*cos(omega(n)*t) + phidotg(n)*sin(omega(n)*t)/omega(n);
      Psi(k,i) = Psi(k,i) + phi(n,i) * buf;
    end
  end
end

				% GRAPH OF EIGENFUNCTIONS

if (strcmp(graph_eigenfunctions,'yes')==1); 
  for n=1:nmax
    g1=plot(x/L,phi(n,:));
    hold on;
  end
  xlim([xmin/L,xmax/L]);
  ylim([-sqrt(2/L),sqrt(2/L)]);
  line([xmin/L,xmax/L],[0 0],'linestyle','-','color',[0.5 0.5 0.5]); % y=0 line in grey
  waitfor(g1);  
end

				% GRAPH OF INITIAL CONDITIONS

if (strcmp(graph_pinch,'yes')==1);
  g2=plot(x/L,f);
  hold on;
  g2=plot(x/L,g); 
  hold on;  
  xlim([xmin/L,xmax/L]);
  ylim([-z,z]);
  xlabel('x/L');
  ylabel('Amplitude');
  line([xmin/L,xmax/L],[0 0],'linestyle','-','color',[0.5 0.5 0.5]); 
  waitfor(g2);  
end

% GRAPH TO APPRECIATE THE APPROXIMATION OF THE EIGENFUNCTION EXPANSION

if (strcmp(graph_approx,'yes')==1);
  g3=plot(x/L,f,'linestyle',':');
  hold on;
  g3=plot(x/L,fapprox,'linestyle','-');
  hold on;
  
  g3=plot(x/L,g,'linestyle',':');
  hold on;
  g3=plot(x/L,gapprox,'linestyle','-');
  
  xlim([xmin/L,xmax/L]);
  ylim([-z,z]);
  line([xmin/L,xmax/L],[0 0],'linestyle','-','color',[0.5 0.5 0.5]); 
  waitfor(g3);  
end

      % GRAPH OF COMPLETE SOLUTION AT SOME SAMPLE TIMES WITHIN A PERIOD

if (strcmp(graph_sample_times,'yes')==1);
  for k=1:nst
    g4=plot(x/L,Psi(k,:));
    ts{k}=sprintf('t =%#10.2E s',tmin+(k-1)*timestep);
			     % String array elements between braces {}
    hold on;
  end
  lgd=legend(ts,'location','eastoutside'); % Add legend of each curve
  set(lgd,'xlabel',' ','ylabel',' ');
  
  xlim([xmin/L,xmax/L]);
  ylim([-z,z]);
  waitfor(g4);  
end

				% MOVIE

if (strcmp(graph_movie,'yes')==1);
  normx=x/L;
  tmin=0;
  tmax= 2 * Longest_period;
  
  timestep=(tmax-tmin)/(movie-1);
  
  for i=1:npt
    psi(i)=0;
    for n=1:nmax
      psi(i) = psi(i) + phi(n,i) * phidotf(n);
    end
  end
				% Initialize movie
  
  g5=plot(normx,psi,'XDataSource','normx','YDataSource','psi');
  
  xlim([xmin/L,xmax/L]);
  ylim([-z,z]);
  line([xmin/L,xmax/L],[0 0],'linestyle','-','color',[0.5 0.5 0.5]); 

				% Loop on frames

  for k=1:movie
    pause(0); % You may take your time to observe each frame by setting some non zero delay [s]
    t = tmin + (k-1) * timestep;
    for i=1:npt
      psi(i)=0;
      for n=1:nmax
	buf = phidotf(n)*cos(omega(n)*t) + phidotg(n)*sin(omega(n)*t)/omega(n);
	psi(i) = psi(i) + phi(n,i) * buf;
      end
    end
    refreshdata();
  end
 
  waitfor(g5); 
end

