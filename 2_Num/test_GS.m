# 4. Define a squared matrix from 1 to 9
aa_1 = [1, 4, 6; -1, 7, 2; 4, 1, 9]

[m, n] = size(aa_1);

q_1 = zeros(m, n);
r_1 = zeros(n, n);

q_1, r_1 = GS_classic(aa_1, 1);
q_1
r_1


% # 5. Run the script
% sigma = 1E-10;
% matrix_2 = [1, 1, 1; sigma, 0, 0; 0, sigma, 0; 0, 0, sigma]

% out_2 = GS_classic(matrix_2, 1)

% # 6. Hilbert matrices

% for n = 3:10
%     matrix_3 = hilb(n)
%     out_3 = GS_classic(matrix_3, 1)
% end

# 7. Modifications that should be done to GS for complex on C^n, instead of reals on R^n