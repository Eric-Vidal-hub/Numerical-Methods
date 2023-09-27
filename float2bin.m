function b=float2bin(f)

%% CONVERSION OF FLOATING-POINT NUMBER TO BINARY STRING

% Input:  f = single or double precision floating-point number 
% Output: b = string of bits in IEEE 754 floating-point format

% Floating-point binary formats:
% Single precision: 1 sign bit,  8 exponent bits, 23 significand bits
% Double precision: 1 sign bit, 11 exponent bits, 52 significand bits

% Exponent in offset binary-1 representation
% Significand in complement-to-two representation

if (isfloat(f)==0)
  fprintf('From float2bin: Argument is not a floating-point number.');
  return;
end
hex = '0123456789abcdef'; % Hexadecimal characters
h = num2hex(f);	          % Convert float to hexadecimal characters
hc = num2cell(h);         % Convert to cell array of characters
nums = cellfun(@(x) find(hex == x) - 1, hc); % Convert to array of numbers
bins = dec2bin(nums, 4);  % Convert to array of binary number strings
b = reshape(bins.', 1, numel(bins)); % Reshape as horizontal vector

end % End of function b=float2bin(f)