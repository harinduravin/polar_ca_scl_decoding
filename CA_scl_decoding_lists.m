clear

n = 2048;
k = 1024;

r = k/n;

crc_size = 8;
crc_poly = [1,0,1,0,1,0,1,1,1];

load('reliability','reliability_order')

[G_N,G,sys_G,H] = gen_polar_g(n,k,reliability_order);

% Define parameters
list_sizes = [32]; % Different values of L to test
snr_dB = 1:0.25:3;         % SNR range in dB
num_trials = 1e4;         % Number of trials per SNR point

% Frozen bits setup (example)
frozen_bits = zeros(1, n);
frozen_bits(reliability_order(k+1:end)+1) = 1; % Example frozen bits pattern

% Initialize BER array for different list sizes
bler_succ = zeros(length(list_sizes), length(snr_dB));

tic
% Loop over each list size
for l_idx = 1:length(list_sizes)
    L = list_sizes(l_idx); % Current list size
    
    % Loop over each SNR value
    for idx = 1:length(snr_dB)
        snr_dB(idx);
        % Calculate SNR in linear scale
        snr_linear = 10^(snr_dB(idx) / 10);
        noise_variance = 1 / (2 * snr_linear * r); % Noise variance for BPSK in AWGN

        error_count = 0; % Initialize error counter
        total_bits = 0;  % Initialize total bits sent
        bler_count = 0;  % Initialize BLER counter

        % Simulation trials
        for trial = 1:num_trials
            % Generate random message
            info_bits =  randi([0 1], 1, k-crc_size);
            % Append CRC
            % Pad data with zeros for division
            data_padded = [zeros(1, crc_size), info_bits];
            [~, remainder] = gfdeconv(data_padded, crc_poly, 2); % Binary division
            
            % Compute CRC
            if length(remainder) < crc_size
                crc_bits = [remainder, zeros(1, crc_size - length(remainder))];
            else
                crc_bits = remainder;
            end
            
            % Add CRC to the data
            message = [crc_bits,info_bits];
    
            % Encode message using the generator matrix
            codeword = mod(message * G, 2);
    
            % BPSK modulation (0 -> 1, 1 -> -1)
            bpsk_signal = 1 - 2 * codeword; 
    
            % Add AWGN noise
            noise = sqrt(noise_variance) * randn(1, n);
            received_signal = bpsk_signal + noise;
           
            [dec_msg] = SCLDecode(L, noise_variance, n, frozen_bits, received_signal, G_N, crc_poly);

            % Count errors
            error_count = error_count + sum(dec_msg ~= message);

            if any(dec_msg ~= message)
                bler_count = bler_count + 1;
            end
            total_bits = total_bits + k;
        end

        % Calculate BLER
        bler_succ(l_idx, idx) = bler_count / num_trials;
        fprintf('List size: %d, SNR (dB): %d, BLER: %g\n', L, snr_dB(idx), bler_succ(l_idx, idx));
    end
end
toc
% Plot BLER vs SNR for different list sizes

% figure;
% for l_idx = 1:length(list_sizes)
%     semilogy(snr_dB, bler_succ(l_idx, :), '-o', 'LineWidth', 2 , 'DisplayName', sprintf('L = %d', list_sizes(l_idx)))
%     hold on;
% end
% grid on;
% xlabel('SNR (dB)');
% ylabel('Block Error Rate (BLER)');
% title('BLER vs SNR for BPSK with SCL Decoder in AWGN Channel');
% legend show;
% hold off;

save('ca_mul_list_size')

% Upper Level Functions

function [decoded_msg] = SCLDecode(L, noise_var, n, frozen_bits, y, G_N, crc_poly)
    m = log2(n); % m is the log2 of the code length

    [inactivePathIndices, activePath, arrayPointer_LLR, arrayPointer_C, ...
        pathIndexToArrayIndex, inactiveArrayIndices, arrayReferenceCount, LLR_path_metric] = initializeDataStructures(L, m);
    [l,inactiveArrayIndices, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount] = ...
    assignInitialPath(inactiveArrayIndices, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, m);
    [~, LLR_0, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
        getArrayPointer_LLR(0, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);


    % Initialization
    for beta = 0:n-1
        % P(1,beta+1, 1) = W(1, y(beta+1)+1); % Assuming W is a matrix with likelihoods
        % P(1,beta+1, 2) = W(2, y(beta+1)+1);
        LLR_0(beta+1) = 2*y(beta+1)/noise_var;
    end

    [arrayPointer_LLR] = setArrayPointer_LLR(0, l, LLR_0 , arrayPointer_LLR, pathIndexToArrayIndex);

    % Main decoding loop
    for phi = 0:n-1

        [inactiveArrayIndices, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex]...
            = recursivelyCalcLLR(m, phi, m, L, inactiveArrayIndices, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, activePath); % Calculate LLRs

        if frozen_bits(phi+1) % Check if u_hat is frozen
            for l = 0:L-1
                if ~activePath(l + 1) % MATLAB indices are 1-based
                    continue;
                end
                [S_p, C_m, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
                    getArrayPointer_C(m, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
                C_m(1,mod(phi,2)+1) = 0; % Set to frozen value
                [arrayPointer_C] = setArrayPointer_C(m, l, C_m , arrayPointer_C, pathIndexToArrayIndex);
                [~,LLR_m, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
                    getArrayPointer_LLR(m, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
                LLR_path_metric(S_p+1) = LLR_path_metric(S_p+1) + log(1+exp(-LLR_m(1)));
            end
            
        else
            [activePath, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, inactiveArrayIndices, inactivePathIndices, LLR_path_metric]...
            = continuePaths_UnfrozenBit(phi, L, m, activePath, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, inactiveArrayIndices, inactivePathIndices, LLR_path_metric);
        end

        if mod(phi, 2) == 1
            [inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex]...
                = recursivelyUpdateC(m, phi, m, L, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex, activePath);
        end
    end

    l_prime = -1;
    p_prime = realmax;

    crc_passed = 0;

    for l = 0:L-1
        [~,C_0, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
            getArrayPointer_C(0, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);

        est_codeword = C_0(:,1).';

        check_u = mod(est_codeword*G_N,2);
        check_msg = check_u(not(frozen_bits));

        [~, check_remainder] = gfdeconv(check_msg, crc_poly, 2);
        
        % Check if remainder is all zeros
        if ~all(check_remainder == 0)
            continue;
        end

        % disp(l)
        crc_passed = 1;
        if p_prime > LLR_path_metric(l+1)
            l_prime = l;
            p_prime = LLR_path_metric(l+1);
        end

    end

    % fprintf('done crc.\n');
    if crc_passed == 0
        % fprintf('Starting normal one.\n');
        for l = 0:L-1
            if ~activePath(l + 1) % MATLAB indices are 1-based
                continue;
            end

            if p_prime > LLR_path_metric(l+1)
                l_prime = l;
                p_prime = LLR_path_metric(l+1);
            end
        end
        % disp(l_prime)
        % fprintf('done normal.\n');
    end

    [~, C_0, ~, ~, ~, ~, ~] = ...
          getArrayPointer_C(0, l_prime, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
    
    est_codeword = C_0(:,1).'; % Return the decoded word
    check_u = mod(est_codeword*G_N,2);
    decoded_msg = check_u(not(frozen_bits));

    % disp(l_prime)
    % fprintf('done final.\n');
end

function [activePath, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, inactiveArrayIndices, inactivePathIndices, LLR_path_metric]...
= continuePaths_UnfrozenBit(phi, L, m, activePath, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, inactiveArrayIndices, inactivePathIndices, LLR_path_metric)
    % Initialize probForks and its iterator
    LLRForks = -inf(L, 2); % Initialize with 0
    i = 0;

    % Populate LLRForks
    for l = 0:(L - 1)
        if activePath(l + 1) % Check if path `l` is active
            [S_p, L_m, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
                getArrayPointer_LLR(m, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
            LLRForks(l+1,1) = -(LLR_path_metric(S_p+1) + log(1+exp(-L_m(1))));
            LLRForks(l+1,2) = -(LLR_path_metric(S_p+1) + log(1+exp(L_m(1))));
            i = i + 1;
        end
    end

    % Determine threshold rho
    rho = min(2 * i, L);

    % Initialize contForks
    % Populate contForks such that the top `rho` probabilities are selected
    allForks = reshape(LLRForks, [], 1); % Flatten into 1D array
    [~, indices] = maxk(allForks, rho); % Find indices of top `rho` elements
    contForks = zeros(size(allForks));
    contForks(indices) = true;
    contForks = reshape(contForks,size(LLRForks));

    % First, kill-off non-continuing paths
    for l = 0:(L - 1)
        if ~activePath(l + 1) % Skip inactive paths
            continue;
        end
        if ~contForks(l + 1, 1) && ~contForks(l + 1, 2)
            [inactivePathIndices, activePath, inactiveArrayIndices, arrayReferenceCount, LLR_path_metric] = ...
            killPath(l, inactivePathIndices, activePath, pathIndexToArrayIndex, inactiveArrayIndices, arrayReferenceCount, LLR_path_metric, m); 
            % Kill path if both forks are bad
        end
    end

    % Then, continue relevant paths and duplicate if necessary
    for l = 0:(L - 1)
        if ~activePath(l + 1) % Skip inactive paths
            continue;
        end
        if ~contForks(l + 1, 1) && ~contForks(l + 1, 2)
            continue; % Skip bad forks
        end

        % Get the pointer to bit pairs
        [~,Cm, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
            getArrayPointer_C(m, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex); % Get the pointer to bit pairs
        [~, L_m_l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
            getArrayPointer_LLR(m, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);

        if contForks(l + 1, 1) && contForks(l + 1, 2)
            % Both forks are good
            Cm(1,mod(phi,2)+1) = 0;

            % Setter
            [arrayPointer_C] = setArrayPointer_C(m, l, Cm , arrayPointer_C, pathIndexToArrayIndex);

            [l_prime, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, LLR_path_metric] = ...
                clonePath(l, m, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, LLR_path_metric);

            [~,Cm, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
                getArrayPointer_C(m, l_prime, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
            Cm(1,mod(phi,2)+1) = 1; % Set for cloned path

            % Setter
            [arrayPointer_C] = setArrayPointer_C(m, l_prime, Cm , arrayPointer_C, pathIndexToArrayIndex);

            [~, L_m_l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
                getArrayPointer_LLR(m, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
            [~, L_m_l_prime, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
                getArrayPointer_LLR(m, l_prime, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
            LLR_path_metric(l+1) = LLR_path_metric(l+1) + log(1+exp(-L_m_l(1)));
            LLR_path_metric(l_prime+1) = LLR_path_metric(l_prime+1) + log(1+exp(L_m_l_prime(1)));

        elseif contForks(l + 1, 1)
            % Only fork 0 is good
            Cm(1,mod(phi,2)+1) = 0;
            [arrayPointer_C] = setArrayPointer_C(m, l, Cm , arrayPointer_C, pathIndexToArrayIndex);
            LLR_path_metric(l+1) = LLR_path_metric(l+1) + log(1+exp(-L_m_l(1)));
        elseif contForks(l + 1, 2)
            % Only fork 1 is good
            Cm(1,mod(phi,2)+1) = 1;
            [arrayPointer_C] = setArrayPointer_C(m, l, Cm , arrayPointer_C, pathIndexToArrayIndex);
            LLR_path_metric(l+1) = LLR_path_metric(l+1) + log(1+exp(L_m_l(1)));
        end

    end
end

% Mid level functions

% List version
function [inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex]...
    = recursivelyUpdateC(lambda, phi, m, L, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex, activePath)
    if mod(phi, 2) == 0
        return;
    end
    
    psi = floor(phi / 2);
    for l = 0:L-1
        if ~activePath(l + 1) % MATLAB indices are 1-based
            continue;
        end
        [~,C_lamb, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = getArrayPointer_C(lambda, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
        [~,C_lamb_1, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = getArrayPointer_C(lambda-1, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
        for beta = 0:(2^(m-lambda))-1
    
            % Update C values using the XOR operation
            C_lamb_1(2*beta+1,mod(psi,2)+1) = mod(C_lamb(beta+1,1)+ C_lamb(beta+1,2),2);
            C_lamb_1(2*beta+1+1,mod(psi,2)+1) = C_lamb(beta+1,2);
        end
        % Setters
        [arrayPointer_C] = setArrayPointer_C(lambda-1, l, C_lamb_1 , arrayPointer_C, pathIndexToArrayIndex);
    end

    if mod(psi, 2) == 1
        [inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex]...
            = recursivelyUpdateC(lambda-1, psi, m, L, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex, activePath);
    end
end

% List version
function [inactiveArrayIndices, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex]...
    = recursivelyCalcLLR(lambda, phi, m, L, inactiveArrayIndices, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, activePath)
    if lambda == 0
        return; % Base case
    end

    psi = floor(phi / 2);
    if mod(phi, 2) == 0
        [inactiveArrayIndices, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex]...
            = recursivelyCalcLLR(lambda - 1, psi, m, L, inactiveArrayIndices, arrayPointer_LLR, arrayPointer_C, arrayReferenceCount, pathIndexToArrayIndex, activePath);
    end

    sigma = 0;
    for l = 0:L-1
        if ~activePath(l + 1) % MATLAB indices are 1-based
            continue;
        end
        [~,LLR_lamb, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = getArrayPointer_LLR(lambda, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
        [~,LLR_lamb_1, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = getArrayPointer_LLR(lambda-1, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);
        [~,C_lamb, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = getArrayPointer_C(lambda, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex);

        for beta = 0:(2^(m-lambda))-1
            if mod(phi, 2) == 0
                % Apply Equation (1)
                LLR_lamb(beta+1) = cnop_llr(LLR_lamb_1(2*beta+1),LLR_lamb_1(2*beta+1+1));                
            else
                % Apply Equation (2)
                u_prime = C_lamb(beta+1, 1);
                LLR_lamb(beta+1) = vnop_llr((1-2*u_prime)*LLR_lamb_1(2*beta+1),LLR_lamb_1(2*beta+1+1));                
            end
        end

        % Setter
        [arrayPointer_LLR] = setArrayPointer_LLR(lambda, l, LLR_lamb , arrayPointer_LLR, pathIndexToArrayIndex);
    end

end


% Low level functions
function [inactivePathIndices, activePath, arrayPointer_LLR, arrayPointer_C, ...
    pathIndexToArrayIndex, inactiveArrayIndices, arrayReferenceCount, LLR_path_metric] = initializeDataStructures(L, m)

    inactivePathIndices = []; % Stack with capacity L
    activePath = false(1, L); % Boolean array of size L
    
    % 2D arrays of size (m + 1) x L
    arrayPointer_LLR = cell(m + 1, L);
    LLR_path_metric = zeros(L,1);
    arrayPointer_C = cell(m + 1, L);
    pathIndexToArrayIndex = zeros(m + 1, L);
    
    % Array with stacks
    inactiveArrayIndices = cell(1, m + 1);
    arrayReferenceCount = zeros(m + 1, L);
    
    % Initialize arrays
    for lambda = 0:m
        inactiveArrayIndices{lambda + 1} = 0:L-1;
        for s = 1:L
            arrayPointer_LLR{lambda + 1, s} = nan(2^(m - lambda),1); % float pairs
            arrayPointer_C{lambda + 1, s} = nan(2^(m - lambda),2); % bit pairs
            arrayReferenceCount(lambda + 1, s) = 0;
        end
    end
    
    % Initialize active paths
    for l = 0:L-1
        activePath(l+1) = false;
        inactivePathIndices(end + 1) = l;
    end
end

function [l, inactiveArrayIndices, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount] = ...
    assignInitialPath(inactiveArrayIndices, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, m)

    l = inactivePathIndices(end);
    inactivePathIndices(end) = []; % Pop
    activePath(l+1) = true;

    for lambda = 0:m
        s = inactiveArrayIndices{lambda + 1}(end);
        inactiveArrayIndices{lambda + 1}(end) = []; % Pop
        
        pathIndexToArrayIndex(lambda + 1, l+1) = s; 
        arrayReferenceCount(lambda + 1, s+1) = 1;
    end
end

function [l_prime, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, LLR_path_metric] = ...
    clonePath(l, m, inactivePathIndices, activePath, pathIndexToArrayIndex, arrayReferenceCount, LLR_path_metric)

    l_prime = inactivePathIndices(end);
    inactivePathIndices(end) = [];
    activePath(l_prime+1) = true;
    LLR_path_metric(l_prime+1) = LLR_path_metric(l+1);
    
    for lambda = 0:m
        s = pathIndexToArrayIndex(lambda + 1, l+1);
        pathIndexToArrayIndex(lambda + 1, l_prime+1) = s;
        arrayReferenceCount(lambda + 1, s+1) = arrayReferenceCount(lambda + 1, s+1) + 1;
    end
end


function [inactivePathIndices, activePath, inactiveArrayIndices, arrayReferenceCount, LLR_path_metric] = ...
    killPath(l, inactivePathIndices, activePath, pathIndexToArrayIndex, inactiveArrayIndices, arrayReferenceCount, LLR_path_metric, m)

    activePath(l+1) = false;
    inactivePathIndices(end + 1) = l;

    LLR_path_metric(l+1) = 0;

    for lambda = 0:m
        s = pathIndexToArrayIndex(lambda + 1, l+1);
        arrayReferenceCount(lambda + 1, s+1) = arrayReferenceCount(lambda + 1, s+1) - 1;

        if arrayReferenceCount(lambda + 1, s+1) == 0
            inactiveArrayIndices{lambda + 1}(end + 1) = s; % Push
        end
    end
end

function [s_p, pointer, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
    getArrayPointer_LLR(lambda, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex)

    s = pathIndexToArrayIndex(lambda+1, l+1);
    if arrayReferenceCount(lambda+1, s+1) == 1
        pointer = arrayPointer_LLR{lambda + 1, s+1};
        s_p = s;
    else
        s_prime = inactiveArrayIndices{lambda + 1}(end);
        inactiveArrayIndices{lambda + 1}(end) = [];
        
        arrayPointer_LLR{lambda + 1, s_prime+1} = arrayPointer_LLR{lambda + 1, s+1}; % Copy contents
        arrayPointer_C{lambda + 1, s_prime+1} = arrayPointer_C{lambda + 1, s+1}; % Copy contents
        arrayReferenceCount(lambda + 1, s+1) = arrayReferenceCount(lambda + 1, s+1) - 1;
        arrayReferenceCount(lambda + 1, s_prime+1) = 1;
        pathIndexToArrayIndex(lambda + 1, l+1) = s_prime;
        
        pointer = arrayPointer_LLR{lambda + 1, s_prime + 1};
        s_p = s_prime;
    end
end

function [s_p, pointer, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex] = ...
    getArrayPointer_C(lambda, l, inactiveArrayIndices, arrayPointer_C, arrayPointer_LLR, arrayReferenceCount, pathIndexToArrayIndex)

    s = pathIndexToArrayIndex(lambda+1, l+1);

    if arrayReferenceCount(lambda+1, s+1) == 1
        pointer = arrayPointer_C{lambda + 1, s+1};
        s_p = s;
    else
        s_prime = inactiveArrayIndices{lambda + 1}(end);
        inactiveArrayIndices{lambda + 1}(end) = [];
        
        arrayPointer_C{lambda + 1, s_prime+1} = arrayPointer_C{lambda + 1, s+1}; % Copy contents
        arrayPointer_LLR{lambda + 1, s_prime+1} = arrayPointer_LLR{lambda + 1, s+1}; % Copy contents
        arrayReferenceCount(lambda + 1, s+1) = arrayReferenceCount(lambda + 1, s+1) - 1;
        arrayReferenceCount(lambda + 1, s_prime+1) = 1;
        pathIndexToArrayIndex(lambda + 1, l+1) = s_prime;
        
        pointer = arrayPointer_C{lambda + 1, s_prime + 1};
        s_p = s_prime;
    end
end

function [arrayPointer_LLR] = setArrayPointer_LLR(lambda, l, pointer , arrayPointer_LLR, pathIndexToArrayIndex)

    s = pathIndexToArrayIndex(lambda+1, l+1);
    arrayPointer_LLR{lambda + 1, s + 1} = pointer;
end

function [arrayPointer_C] = setArrayPointer_C(lambda, l, pointer , arrayPointer_C, pathIndexToArrayIndex)

    s = pathIndexToArrayIndex(lambda+1, l+1);
    arrayPointer_C{lambda + 1, s + 1} = pointer;
end

function l = cnop_llr(l1,l2)
    l = 2 * atanh(tanh(l1/2).*tanh(l2/2));
end

function l = vnop_llr(l1,l2)
    l = l1 + l2;
end

