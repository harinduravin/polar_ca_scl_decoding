
function [G,gen_matrix,sys_gen_matrix, parity_matrix] = gen_polar_g(N,K,reliability_order)

% Parameters
% reliability_order = [3,5,6,7,8,1,2,4]-1;
% reliability_order =  7:-1:0;
% reliability_order = [127 126 125 123 119 111 95 63 124 122 121 118 117 110 115 109 107 94 103 93 91 87 62 61 79 59 55 47 31 120 116 114 108 113 106 105 102 92 101 90 89 99 86 85 60 83 78 58 77 57 54 75 53 51 46 71 45 43 30 39 29 27 112 104 100 23 88 98 84 97 82 76 56 81 15 74 52 73 50 70 44 49 69 42 41 67 38 28 37 26 25 35 22 96 21 80 14 19 72 13 48 68 40 11 66 36 65 24 7 34 20 33 18 12 17 10 64 9 6 5 32 3 16 8 4 2 1 0];
% Reliability order from most to least reliable

% Step 1: Define frozen and information bit indices
frozen_indices = sort(reliability_order(K+1:end)); % Less reliable positions (freeze)
info_indices = sort(reliability_order(1:K)); % Most reliable positions (information bits)

% Step 2: Initialize the base polar transformation matrix (N x N)
F = [1 0; 1 1]; % Basic 2x2 kernel
F_N = 1; % Start with identity

% Kronecker product to build the full matrix
for i = 1:log2(N)
    F_N = kron(F_N, F);
end

B_N = generate_polar_permutation_matrix(N);

G = mod(B_N*F_N,2);

% Step 3: Construct the generator matrix for Polar code
% Only keep the rows of G corresponding to information bit positions
gen_matrix = G(info_indices + 1, :); % +1 to account for MATLAB 1-based indexing
parity_matrix = G(:, frozen_indices + 1);
parity_matrix = parity_matrix.';

[permuted_info_indices, comp_permuted_info_indices] = permute_indices(B_N, info_indices + 1);

G_AB = G(info_indices + 1, permuted_info_indices);
G_ABc = G(info_indices + 1, comp_permuted_info_indices);

P_G = inv(gf(G_AB))*gf(G_ABc);

G_sys = [eye(K) P_G];

sys_gen_matrix = double(G_sys.x);


% Display the generator matrix
% disp('Generator matrix for Polar code:')
% disp(G_polar)

end

function B_N = generate_polar_permutation_matrix(N)
    % Generate the permutation matrix B_N for polar codes of size N.
    % N must be a power of 2 (e.g., 2, 4, 8, 16, ...).

    if mod(log2(N), 1) ~= 0
        error('N must be a power of 2.');
    end

    % Start with the base case B_1, which is simply [1]
    B_N = 1;

    % Recursive generation of B_N using the formula
    for i = 1:log2(N)
        R_N = reverse_order(2^i);
        B_N = mod(R_N * kron(eye(2), B_N),2);
    end
end

function R = reverse_order(N)
    % Generate the reverse ordering matrix R_N of size N
    R = zeros(N);
    for i = 1:N
        % Bit-reversal index for position i-1
        reversed_index = bitrevorder(i-1, N) + 1;
        R(i, reversed_index) = 1;
    end
end

function r = bitrevorder(index, N)
    % Perform bit-reversal of index with log2(N) bits
    num_bits = log2(N);
    r = 0;
    % Shift the bits left by one position, wrapping the leftmost bit to the rightmost position
    for i = 0:num_bits-1
        if bitand(index, bitshift(1, i))
            % Rotate left by one: bit `i` maps to `i-1` (wrapping around)
            new_pos = mod(i-1, num_bits);
            r = bitor(r, bitshift(1, new_pos));
        end
    end
end

function [permuted_indices,comp_permuted_indices] = permute_indices(P, indices)
    % Obtain the image of a set of row indices through a permutation matrix
    % Inputs:
    %   P - An N x N permutation matrix
    %   indices - A vector of row indices (1-based) to be permuted
    % Output:
    %   permuted_indices - The indices of the rows after applying P

    % Check if P is a valid permutation matrix
    if ~isequal(P * P', eye(size(P))) || ~isequal(P' * P, eye(size(P)))
        error('P must be a permutation matrix.');
    end
    
    % Convert indices into a one-hot representation to apply P
    N = size(P, 1);
    one_hot_indices = zeros(N, 1);
    one_hot_indices(indices) = 1;
    
    % Apply the permutation matrix to the one-hot vector
    permuted_one_hot = P * one_hot_indices;
    
    % Find the indices that are set to 1 in the permuted one-hot vector
    permuted_indices = find(permuted_one_hot == 1);
    comp_permuted_indices = find((1-permuted_one_hot) == 1);
end


