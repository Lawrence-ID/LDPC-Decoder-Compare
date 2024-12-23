function [intera, EbN0, BER_flooding, BER_layered, FER_flooding, FER_layered] = read_statistic(filename, LDPC_K)
    iteras_str = regexp(filename, '(\d+)iteras', 'match');
    intera_str = iteras_str{1}(1:end-6); % 去除 'iteras' 后缀
    intera = str2double(intera_str);
    result = importdata(filename);

    EbN0 = result.data(:, 1);
    BER_flooding = result.data(:, 2);
    BER_layered = result.data(:, 3);
    FER_flooding = 1 - (1 - BER_flooding) .^ LDPC_K;
    FER_layered = 1 - (1 - BER_layered) .^ LDPC_K;
end