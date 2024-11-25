clear

% Parameters
N = 2048;        % Block length (must be a power of 2)
K = 1024;         % Number of information bits
SNR_dB = 3;   % Signal-to-noise ratio in dB

% Convert SNR from dB to linear scale
SNR = 10^(SNR_dB / 10);

% SNR = 1/1.261914;

channels = zeros(1,N);
block_err = zeros(1,N);

for i = 1:N
    channels(i) = get_exp_llr(SNR,N,i);
    block_err(i) = 0.5*erfc(0.5*sqrt(channels(i)));
end

% Plot with star markers and no connecting lines
semilogy(1:N, block_err(1:N), '*', 'MarkerSize', 6);
xlabel('Channel Index');
ylabel('Block Error Probability');
title('Block Error Probability for Polar Code Channels');
grid on;

[~, sortedIndices] = sort(block_err);
reliability_order = sortedIndices-1;
save('reliability','reliability_order')

% uniqueValues = setxor(sortedIndices(1:116)-1, reliability_order(1:116));
% uniqueValues

function llr = get_exp_llr(SNR,N,i)

    if N == 1 && i == 1
        llr = 2*SNR;
    elseif N == 1 && i~= 1
        disp('cat')
    elseif mod(i,2) == 0
        llr = 2*get_exp_llr(SNR,N/2,i/2);
    else
        llr = phi_inv(1 - (1 - phi(get_exp_llr(SNR,N/2,(i+1)/2)))^2);        
    end
end


% Initialize the reliability of channels
% z = zeros(1, N);  
% z(1) = 2 * SNR;   

% Calculate the reliability of each bit channel
% for n = 1:log2(N)
%     % Loop through pairs of channels
%     for j = 1:2^(n-1)
%         T = z(j);
%         z(j) = phi_inv(1 - (1 - phi(T))^2);   
%         z(2^(n-1) + j) = 2 * T;               
%     end
% end


function phi_val = phi(x)
    if x < 10
        phi_val = exp(-0.4527 * x^0.86 + 0.0218);
    else
        phi_val = sqrt(pi/x)*exp(-x/4)*(1-10/(7*x));
    end
end

function x = phi_inv(y)
    x = (-(log(y) - 0.0218) / 0.4527)^(1/0.86);
end
