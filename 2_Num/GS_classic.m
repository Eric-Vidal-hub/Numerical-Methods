function [q, r] = GS(aa, normalize)

% get the size of the input matrix
[m, n] = size(aa);

if m<n
    fprintf('The number of rows must be greater or equal to the number of columns');
    stop();
end

% initialize the q and r matrices
q = zeros(m, n);
r = zeros(m, n);

for j=1:n
    q(:,j) = aa(:, j);
    for i=1:j-1
        r(i,j) = q(:,i)' * aa(:,j);
        q(:,j) = q(:,j) - r(i,j) * q(:,i);
    end
end