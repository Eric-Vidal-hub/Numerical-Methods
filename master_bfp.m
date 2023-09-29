clear all;  % Clear all variables
clc;        % Clear the command window

% 3.4.1 Mastering floating-point numbers

% 1. Evaluation of the numbero fo significant decimal digits for normalized bfp numbers
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
    bit_pattern = dec2bin(typecast(u_k, 'uint64'));
    fprintf('%2d %15.7E %4d %s\n', k, u_k, dec_dig_k, bit_pattern)
end
