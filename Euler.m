#! /opt/local/bin/octave
clear all;   % Clear all variables/functions in memory

% Euler.m

pi           % pi is predefined in Matlab/Octave
ci = 1i      % Imaginary number

graph='yes'; % Character variable

fmt=3;       % Values to test = 1,2,3,4

if (fmt==1); % Testing numerical variables
  fmt_string = '%d %d %d \n';
elseif (fmt==2);
  fmt_string = '%15.6E %15.6E %15.6E \n';
elseif (fmt==3);
  fmt_string = '%#15.6G %#15.6G %#15.6G \n';
end

n = 101 ;    %  Number of data points
xmin = 0;    %  Minimum abscissa
xmax = 2*pi; %  Maximum abscissa
step = (xmax-xmin)/(n-1);

for k=1:n
  x(k) = xmin + (k-1) * step ; % Compute abscissas
end

z = exp(ci*x); % Compute ordinates

f1 = fopen ('Euler.dat','w'); % Open output file
for k=1:n
  fprintf(fmt_string,x(k),real(z(k)),imag(z(k))); % Write to file that you put first which would be f1, if there is no file it is printed at the terminal
end
fclose(f1); % Close output file

if (strcmp(graph,'yes')==1);   % Testing identical character strings
  g=plot(x,real(z),x,imag(z)); % Simple plot in graphical window
  waitfor(g); % Wait for closing graphical window
end

% COMMENTS

% "Octave/Matlab output precision" should not be confused
% with the number of SIGNIFICANT DIGITS.

% To correctly manage significant digits as required in physics and engineering, use:

% (1) SCIENTIFIC NOTATION

% fprintf (' %#w.mG ',y);  % where m   is the number of significant digits

% OR

% (2) NORMALIZED SCIENTIFIC NOTATION (always one non-zero digit before decimal point)

% fprintf (' %w.mE ',y);   % where m+1 is the number of significant digits.
% In this case, m is the "Octave/Matlab output precision",
% i.e. the number of digits after the radix (decimal point).

% In both cases (1) and (2):
% y = real number to be printed respecting significant digit rule.
% w = integer (units of character count) = width of field to print y.
% Typing g instead of G or e instead  E in the format definition toggles
%   the exponent flag respectively in lower or upper case (e or E).

% Both above formats (1) and (2) force to print trailing zeros
% that are significant digits and are useful to perform alignment of columns of numbers.

% MANDATORY : w > m+7 because :
%
% 1 character  for the sign of y
% 1 character  for the decimal point
% 1 character  for the exponent flag (e or E)
% 1 character  for the exponent sign
% 3 characters for the exponent (integer number)

% A good choice is w>m+8 to warrant at least one blank character
% between two numbers on a same line.
