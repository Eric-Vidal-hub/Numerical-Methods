%% Spatial discretisation grids
% generates spatial discretisation grids for problem
directunit = 1/recipunit;   % Unit of direct space

%% Correction
xres = 0.5;     % Ang                 % Resolution of x grid
xmin = -20;     % Ang                 % start of xgrid
xmax = 100;     % Ang                 % end of xgrid
Uxmin1 = 0;     % Ang                 % start of 1st barrier       
Uxmax1 = 15;    % Ang                 % end of 1st barrier
Uxmin2 = 65;    % Ang                 % start of 2nd barrier
Uxmax2 = 80;    % Ang                 % end of 2nd barrier

% shifting all start points by a factor of half resolution to block x=0
% bring included in the grid
xmin = xmin + xres/2;
xmax = xmax + xres/2;
Uxmin1 = Uxmin1 + xres/2;
Uxmin2 = Uxmin2 + xres/2;
Uxmax1 = Uxmax1 - xres/2;
Uxmax2 = Uxmax2 - xres/2;

x = ((xmin:xres:(xmax)))';
L1 = length(x);

Ux = (Uxmin1:xres:(Uxmax2))';    % Potential barrier discretisation
L2 = length(Ux);

% sanity check for limits
lowlim = min([Uxmin1 Uxmin2 Uxmax1 Uxmax2]);
hilim = max([Uxmin1 Uxmin2 Uxmax1 Uxmax2]);

if xmin>lowlim
    fprintf('\n WARNING: Lower limit of x (%d) is greater than start of potential step %d\n',xmin,lowlim)
end

if xmax<hilim
    fprintf('\n WARNING: Upper limit of x (%d) is smaller than end of potential step %d\n',xmax,hilim)
end

fprintf('\n Spatial Discretisation Resolution = %5d\n',xres)