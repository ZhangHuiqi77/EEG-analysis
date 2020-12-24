function [data] = touv(data)
    data = data * 4.5/(2^23 - 1) * 10^6; % transfer the value to uV
end
