clear
info_length = 116;
% n = 128;
crc_size = 8;
crc_poly = [1,0,1,0,1,0,1,1,1];

info_bits =  randi([0 1], 1, info_length-crc_size);

% info_bits = flip([1,1,0,1,0,0,1,1,1,0,1,1,0,0]);


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
transmitted_bits = [crc_bits,info_bits];
fprintf('Transmitted Bits: %s\n', num2str(transmitted_bits));

% Receiver Side
% Extract data and CRC


% Check CRC correctness
% Pad received data with zeros for division
received_padded = [transmitted_bits];
[~, check_remainder] = gfdeconv(received_padded, crc_poly, 2);

% Check if remainder is all zeros
if all(check_remainder == 0)
    fprintf('CRC Check Passed: Data is Correct.\n');
else
    fprintf('CRC Check Failed: Data is Corrupted.\n');
end

