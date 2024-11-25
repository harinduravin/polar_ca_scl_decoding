clear

n = 8;
k = 5;

r = k/n;

reliability_order = 7:-1:0;
% reliability_order = [127 126 125 123 119 111 95 63 124 122 121 118 117 110 115 109 107 94 103 93 91 87 62 61 79 59 55 47 31 120 116 114 108 113 106 105 102 92 101 90 89 99 86 85 60 83 78 58 77 57 54 75 53 51 46 71 45 43 30 39 29 27 112 104 100 23 88 98 84 97 82 76 56 81 15 74 52 73 50 70 44 49 69 42 41 67 38 28 37 26 25 35 22 96 21 80 14 19 72 13 48 68 40 11 66 36 65 24 7 34 20 33 18 12 17 10 64 9 6 5 32 3 16 8 4 2 1 0];

[G_N,G,sys_G,H] = gen_polar_g(n,k);

% sys_H = gen2par(sys_G);


%%%%%%%%%%%%%%%%
% message = randi([0 1], 1, k);
% 
% % Encode message using the generator matrix
% codeword = mod(message * sys_G, 2);
% 
% paritycheck = mod(codeword*(H.'),2);

%%%%%%%%%%%%%%%%

% p = 0.1;

% W =[1-p p;p 1-p]; % Random likelihoods (example)

% Generate random message
% message = randi([0 1], 1, k);
% message = [1 1 1 0 0];

% Encode message using the generator matrix
% codeword = mod(message * G, 2);

% flip_mask = rand(1, n) < p; % Generate a mask where each bit has a probability p to be flipped


% received_codeword = mod(codeword + flip_mask, 2); % Apply bit-flips where flip_mask is 1
% 
frozen_bits = zeros(1,n);
frozen_bits(reliability_order(k+1:end)+1) = 1; % Example frozen bits pattern
% [decoded_word, P, B] = SCDecode(W, n, frozen_bits,received_codeword);
% disp('Error pattern:');
% disp(flip_mask);
% disp('Success:');
% disp(all(codeword == decoded_word));
% disp('Vector u:');
% dec_u = mod(decoded_word*G_N,2);
% disp(dec_u);
% disp('Decoded msg:');
% disp(dec_u(not(frozen_bits)));
% B(3+1,:)


snr_dB = 4:0.5:7;   % SNR range in dB
num_trials = 1e5;  % Number of trials per SNR point

% Initialize BER array
ber_sc = zeros(size(snr_dB));

% Loop over each SNR value
for idx = 1:length(snr_dB)
    snr_dB(idx)
    % Calculate SNR linear scale
    snr_linear = 10^(snr_dB(idx) / 10);
    noise_variance = 1 / (2*snr_linear*r); % Noise variance for BPSK in AWGN

    error_count = 0; % Initialize error counter
    total_bits = 0;  % Initialize total bits sent
    bler_count = 0;

    % Simulation trials
    for trial = 1:num_trials
        % Generate random message
        message = randi([0 1], 1, k);

        % Encode message using the generator matrix
        codeword = mod(message * G, 2);

        % BPSK modulation (0 -> 1, 1 -> -1)
        bpsk_signal = 1 - 2 * codeword; 

        % Add AWGN noise
        noise = sqrt(noise_variance) * randn(1, n);
        received_signal = bpsk_signal+noise;
        
        % Demodulation (hard decision: 0 if > 0, 1 if <= 0)
        % received_word = received_signal < 0;
        
        % Decode using GRAND
        % [decoded_word, success] = GRAND_decoder(sys_G,sys_H, received_word, 1000);
        [decoded_word, L, B] = SCDecode(noise_variance, n, frozen_bits,received_signal);

        dec_u = mod(decoded_word*G_N,2);

        dec_msg = dec_u(not(frozen_bits));
        
        % Count errors
        error_count = error_count + sum(dec_msg ~= message);
        if any(dec_msg ~= message)
            bler_count = bler_count + 1;
        end
        total_bits = total_bits + k;
    end
    
    % Calculate BLER
    bler_succ(idx) = bler_count*k / total_bits;
    fprintf('SNR (dB): %d, BLER: %g\n', snr_dB(idx), bler_succ(idx));
end

% Plot BER curve
figure;
semilogy(snr_dB, bler_succ, '-o', 'LineWidth', 2);
grid on;
xlabel('Eb/No (dB)');
ylabel('Block Error Rate (BLER)');
title('BLER vs SNR for BPSK with SC decoder n = 8 in AWGN Channel');
legend('SC Decoding');


function [decoded_word,L,B] = SCDecode(noise_var, n, frozen_bits,y)
    m = log2(n); % m is the log2 of the code length
    % P = zeros(m+1,n, 2); % Initialize probabilities
    L = zeros(m+1,n); % Initialize LLRs
    B = zeros(m+1,n, 1); % Initialize decoded bits

    % Initialization
    for beta = 0:n-1
        % P(1,beta+1, 1) = W(1, y(beta+1)+1); % Assuming W is a matrix with likelihoods
        % P(1,beta+1, 2) = W(2, y(beta+1)+1);
        L(1,beta+1) = 2*y(beta+1)/noise_var;
    end

    % Main decoding loop
    for phi = 0:n-1
        % P = recursivelyCalcP(m, phi, P, m, B); % Calculate probabilities
        L = recursivelyCalcL(m, phi, L, m, B); % Calculate LLRs

        if frozen_bits(phi+1) % Check if u_hat is frozen
            B(m+1, phi+1) = 0; % Set to frozen value
        else
            % if P(m+1, phi+1, 1) > P(m+1, phi+1, 2)
            %     B(m+1, phi+1) = 0;
            % else
            %     B(m+1, phi+1) = 1;
            % end
            if L(m+1, phi+1) > 0
                B(m+1, phi+1) = 0;
            else
                B(m+1, phi+1) = 1;
            end
        end

        if mod(phi, 2) == 1
            B = recursivelyUpdateB(m, phi, B, m);
        end
    end

    decoded_word = B(1,:); % Return the decoded word
    % dec_ub = B(end,:);
end


function B = recursivelyUpdateB(lambda, phi, B, m)
    if mod(phi, 2) == 0
        return;
    end
    
    psi = floor(phi / 2);
    for beta = 0:(2^(m-lambda))-1
        idx1 = lambdaIndex(lambda-1, psi, 2*beta); 
        idx2 = lambdaIndex(lambda-1, psi, 2*beta+1);

        % Update B values using the XOR operation
        B(lambda-1+1,idx1) = mod(B(lambda+1,lambdaIndex(lambda, phi-1, beta))+ B(lambda+1,lambdaIndex(lambda, phi, beta)),2);
        B(lambda-1+1,idx2) = B(lambda+1,lambdaIndex(lambda, phi, beta));
    end

    if mod(psi, 2) == 1
        B = recursivelyUpdateB(lambda - 1, psi, B, m);
    end
end

function P = recursivelyCalcP(lambda, phi, P, m, B)
    if lambda == 0
        return; % Base case
    end

    psi = floor(phi / 2);
    if mod(phi, 2) == 0
        P = recursivelyCalcP(lambda - 1, psi, P, m, B);
    end

    for beta = 0:(2^(m-lambda))-1
        if mod(phi, 2) == 0
            % Apply Equation (1)
            for u_prime = 0:1

                P(lambda+1, lambdaIndex(lambda, phi, beta), u_prime+1) = ...
                    0.5*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta), mod(u_prime+0,2)+1)*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1), 1)...
                    + 0.5*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta), mod(u_prime+1,2)+1)*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1), 2);

            end
        else
            % Apply Equation (2)
            u_prime = B(lambda+1, lambdaIndex(lambda, phi-1, beta));
            for u_prev = 0:1
                P(lambda+1, lambdaIndex(lambda, phi, beta), u_prev+1) = ...
                    0.5*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta), mod(u_prime+ u_prev,2)+1)*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1), u_prev+1);
            end
        end
    end
end

function L = recursivelyCalcL(lambda, phi, L, m, B)
    if lambda == 0
        return; % Base case
    end

    psi = floor(phi / 2);
    if mod(phi, 2) == 0
        L = recursivelyCalcL(lambda - 1, psi, L, m, B);
    end

    for beta = 0:(2^(m-lambda))-1
        if mod(phi, 2) == 0
            % Apply Equation (1)
            % for u_prime = 0:1
            % 
            %     P(lambda+1, lambdaIndex(lambda, phi, beta), u_prime+1) = ...
            %         0.5*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta), mod(u_prime+0,2)+1)*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1), 1)...
            %         + 0.5*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta), mod(u_prime+1,2)+1)*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1), 2);
            % 
            % end
            L(lambda+1, lambdaIndex(lambda, phi, beta)) = log((custom_exp(L(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta)))*custom_exp(L(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1)))+1)...
                /(custom_exp(L(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta)))+ custom_exp(L(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1)))));

        else
            % Apply Equation (2)
            u_prime = B(lambda+1, lambdaIndex(lambda, phi-1, beta));
            % for u_prev = 0:1
            %     P(lambda+1, lambdaIndex(lambda, phi, beta), u_prev+1) = ...
            %         0.5*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta), mod(u_prime+ u_prev,2)+1)*P(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1), u_prev+1);
            % end
            L(lambda+1, lambdaIndex(lambda, phi, beta)) = log(custom_exp(L(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta)))^(1-2*u_prime)*custom_exp(L(lambda-1+1, lambdaIndex(lambda-1, psi, 2*beta+1))));
        end
    end
end


function idx = lambdaIndex(lambda, phi, beta)
    % Calculate index based on lambda, phi, and beta
    idx = phi + (2^lambda) * beta + 1; % MATLAB index (1-based)
end

function y = custom_exp(x)
    % Set a threshold to avoid underflow
    upper_threshold = 700; % Approximate limit for exp in double precision
    threshold = -700; % Approximate limit for exp in double precision
    y = zeros(size(x)); % Initialize output array
    
    % For values above the threshold, calculate the exponential
    y(x > upper_threshold) = exp(upper_threshold);
    y(x > threshold && x <= upper_threshold) = exp(x(x > threshold && x <= upper_threshold));
    
    % For values below the threshold, use a small approximation to avoid 0
    y(x <= threshold) = exp(threshold);
end


