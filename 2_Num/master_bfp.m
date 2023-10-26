clear all;  % Clear all variables
clc;        % Clear the command window

% 3.4.1 Mastering floating-point numbers

% 1. Evaluation of the number of significant decimal digits for normalized bfp numbers
% DEFINITIONS
% Precision
p(1) = 23;   % Single precision
p(2) = 52;   % Double precision
% Range
q(1) = 8;  % Single precision
q(2) = 11; % Double precision


for i = 1:2
    mant_dig(i) = 2 ^ p(i) - 1;     % Number of significant binary digits, Add phantom digits & use eq. 3.33
    dec_dig(i) = floor(log10(mant_dig(i)));    % Number of significant decimal digits

    % 2. Evaluation of \eps_0 = \eps^{-p} or both sp and dp
    eps0(i) = 2 ^ (-p(i));

    if (i == 1);
        eps0_m(i) = eps(single(1.0)); % Comparison with the built-in function
        
        fprintf('\n SINGLE PRECISION\n');
        fprintf('Number of significant binary digits: %15.7E\n', mant_dig(i));
        fprintf('Number of significant decimal digits: %4d\n', dec_dig(i));
        fprintf('eps_0 = %15.7E\n', eps0(i));
        fprintf('eps_0 (built-in function) = %15.7E\n', eps0_m(i));
    else
        eps0_m(i) = eps(double(1.0)); % Comparison with the built-in function
    
        fprintf('\n DOUBLE PRECISION\n');
        fprintf('Number of significant binary digits: %23.15E\n', mant_dig(i));
        fprintf('Number of significant decimal digits: %4d\n', dec_dig(i));
        fprintf('eps_0 = %23.15E\n', eps0(i));
        fprintf('eps_0 (built-in function) = %23.15E\n', eps0_m(i));
    end
end

% 3. Table creation

% DEFINITIONS
for i = 1:2
    n_min(i) = - 2 ^ (q(i) - 1) + 2;  % Minimum exponent
    n_max(i) = 2 ^ (q(i) - 1) - 1;  % Maximum exponent
    v(i) = (2 - eps0(i)) * 2 ^ n_max(i);     % Number of normalized numbers
    u(i) = 2 ^ n_min(i);                % Underflow gap boundary
    eps_min(i) = 2 ^ (n_min(i) - p(i));       % Minimum relative spacing
end

fprintf('\nTable 1\n')
fprintf('\nPrecision n_min   n_max  eps_0     nu            u        eps_min\n')
fprintf('single:  %5d %5d %12.4E %12.4E %12.4E %12.4E\n', n_min(1), n_max(1), eps0(1), v(1), u(1), eps_min(1))
fprintf('double:  %5d %5d %12.4E %12.4E %12.4E %12.4E\n', n_min(2), n_max(2), eps0(2), v(2), u(2), eps_min(2))

% 4. Intrinsic functions realmin and realmax in dp
fprintf('u-realmin %12.4E\n', u(2) - realmin('double'));
fprintf('v-realmax %12.4E\n', v(2) - realmax('double'));
% Both of them should be 0

% 5. For k till 63, sampling numbers in DP till subnormal

fprintf('\nTable 2\n')
fprintf('\n k      u/2^k   Significant decimal digits    Bit pattern\n')
for k = 1:63
    u_k = u(2) / 2 ^ k;
    dec_dig_k = floor(log10(2 ^ k - 1));
    bit_pattern = float2bin(u_k);
    fprintf('%2d %15.7E %4d %s\n', k, u_k, dec_dig_k, bit_pattern)
end

% % 6. Another table

% fprintf('\nTable 3\n')
% fprintf('\n x(k)=10^k      eps_n(x(k))   eps(x(k))    |x(k)|eps_0    log_10(rho(x(k)))\n')
% for k = -307:308
%     x_k = 10 ^ k;
%     n = floor(log2(x_k));
%     eps_n = 2 ^ (n - p(2));     % Relative spacing for double precision
%     eps_x = eps(x_k);         % Relative spacing for built-in function
%     absx_eps0 = x_k * eps0(2);   % Majoration of epsilon (eq. 3.18)
%     density = -log10(eps_n);    % Log10 Density of floating-point numbers
%     fprintf('%15.7E %15.7E %15.7E %15.7E %15.7E\n', x_k, eps_n, eps_x, absx_eps0, density)
% end

% 7. Another table

fprintf('\nTable 4\n')
fprintf('\n Precision    x=2^{n_{max}+1}   bit pattern\n')
for i = 1:2
    xx(i) = 2 ^ (n_max(i) + 1);
end

fprintf(' Single  %s %s\n', xx(1), float2bin(xx(1)))
fprintf(' Double  %s %s\n', xx(2), float2bin(xx(2)))

% 8. Another table

fprintf('\nTable 5\n')
fprintf('\n Operation    Result   bit pattern\n')
op1 = 0/0
op2 = 0 * Inf
op3 = Inf/Inf
op4 = Inf - Inf
fprintf('0/0 %s %s\n', op1, float2bin(op1))
fprintf('0*Inf %s %s\n', op2, float2bin(op2))
fprintf('Inf/Inf %s %s\n', op3, float2bin(op3))
fprintf('Inf-Inf %s %s\n', op4, float2bin(op4))

% 9. Another table
v = 2 * eps0(2) * 2 ^ n_max(2);  % Max value before overflow in double precision
xm = log(v);                     % Max exponent in an exponential function before overflow in double precision
fprintf('\nTable 6\n')
fprintf('\n x_m      Delta      x = x_m + Delta   exp(x)\n')
for i=0:5
    delta = eps(xm)/2^i;
    x_del = xm + delta;
    fprintf('%15.7E %15.7E %15.7E %15.7E\n', xm, delta, x_del, exp(xm))
end

fprintf('\nTable 7\n')
fprintf('\n x_m      Delta      x = x_m - Delta   exp(x)\n')
for i=-5:1:5
    delta = eps(xm)/2^i;
    x_del = xm - delta;
    fprintf('%15.7E %15.7E %15.7E %15.7E\n', xm, delta, x_del, exp(xm))
end

% 10. Another table

fprintf('\nTable 8\n')
fprintf('\nOffset binary -1\n')
fprintf('\n Precision      n_{min}      n_{max}       u =       v =\n')
for i=1:2
    eps0(i) = 2 ^ (-p(i));
    excess = 2 ^ (q(i) - 1);
    nmax(i) = 2 ^ (q(i)) - 1 - excess - 1;
    nmin(i) = 1 - excess;
    u(i) =
    v(i) = 
    fprintf(' %s %5d %5d %15.7E %15.7E\n', i, n_min(i), n_max(i), u(i), v(i))
end


% 3.4.2 Loss of accuracy in floating-point operations

% Example 1, single precision
x = single(3);
y = single(1E-7);
eps0 = eps(single(1.0));

% Sanity check
if (x<y)
    fprintf('x<y, it doesnt fulfill the requirements.\n')
    stop
end

% Actual code
if ((abs(x) * eps0/2)>y)
    fprintf('Loss of accuracy for single precision in the addition\n')
    fprintf('Results are: %15.7E\n', x+y)
else
    fprintf('No loss of accuracy for single precision in the addition\n')
    fprintf('Results are: %15.7E\n', x+y)
end


% Example 2, double precision
x = 3;
y = 1E-7;
eps0 = eps(1.0);

% Sanity check
if (x<y)
    fprintf('x<y, it doesnt fulfill the requirements.\n')
    stop
end

% Actual code
if ((abs(x) * eps0/2)>y)
    fprintf('Loss of accuracy for double precision in the addition\n')
    fprintf('Results are: %15.7E\n', x+y)
else
    fprintf('No loss of accuracy for double precision in the addition\n')
    fprintf('Results are: %15.7E\n', x+y)
end

% Example 3, single precision

x = 3;
y = 1E-16;
eps0 = eps(1.0);

% Sanity check
if (x<y)
    fprintf('x<y, it doesnt fulfill the requirements.\n')
    stop
end

% Actual code
if ((abs(x) * eps0/2)>y)
    fprintf('Loss of accuracy for double precision in the addition\n')
    fprintf('Results are: %15.7E\n', x+y)
else
    fprintf('No loss of accuracy for double precision in the addition\n')
    fprintf('Results are: %15.7E\n', x+y)
end