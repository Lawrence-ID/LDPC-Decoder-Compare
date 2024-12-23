function [y, llr_final, it] = ldpc_decoder_layered(llr_data, llr_width, intera, HBG, H, check_ctrl, Zc)
    
    llr_ram = llr_data(:).'; % Init llr_ram, convert to row vector

%     HBG=[1 1 1 1 -1 0 1; 0 0 -1 0 0 0 0;2 -1 2 2 2 2 2]; Zc = 3;kb=4;LDPC_K = kb * Zc;LDPC_Q = 3;LDPC_N = Zc * (kb + LDPC_Q);
%     [H, ~, ~] = ldpc_02_HMatrix_calc(HBG, Zc);LDPC_P = Zc * LDPC_Q;  
%     llr_ram = [-10,-34,45,16,-34,30,-60,2,-10,20,-17,-39,12,-5,-14,-7,11,44,4,8,43];


%     HBG=[1 -1 2 1; 0 2 -1 0;2 -1 1 1]; Zc = 3; kb=1; LDPC_K = kb * Zc;LDPC_Q = 3;LDPC_N = Zc * (kb + LDPC_Q);
%     check_ctrl = [3, 3, 3];
%     [H, ~, ~] = ldpc_02_HMatrix_calc(HBG, Zc);
%     llr_ram = [-10,-34,45,16,-34,30,-60,2,-10,20,-17,-39];
%     intera = 2;

    M_c2v_ram = zeros(size(HBG, 1), length(llr_ram));% Init M_c2v_ram: 46 * (68*Zc) * (message bit)

    % HBG: R*C, H: M*N = (R*Zc)*(C*Zc)
    [R_HBG, C_HBG] = size(HBG);
    kb = C_HBG - R_HBG;

    maxMagnitude = 2^(llr_width-1);

    for it = 1 : intera
        for l = 1 : R_HBG

            M_v2c_fifo = zeros(check_ctrl(l), Zc);
            M_v2c_shifted = zeros(check_ctrl(l), Zc);
            read_col_idx = zeros(check_ctrl(l), 1);
            tp = 1;

            % VNUs
            for n = 1 : C_HBG % get Zc LLR at once to VNUs
                if HBG(l, n) ~= -1
                    start_idx = Zc * (n-1) + 1;
                    end_idx = Zc * n;
                    M_v2c_new = llr_ram(start_idx:end_idx) - M_c2v_ram(l, start_idx:end_idx); % M_c2v_l_it-1, 用到上一次迭代的
                    M_v2c_fifo(tp, :) = M_v2c_new;
                    M_v2c_shifted(tp, :) = circshift(M_v2c_new, [0, -HBG(l, n)]); % cyclic shifter(left)
                    read_col_idx(tp, :) = n;
                    tp = tp + 1;
                end
            end

            % 去零、限幅
            M_v2c_shifted(M_v2c_shifted == 0) = 1;%去零
            M_v2c_shifted(M_v2c_shifted > maxMagnitude - 1) = maxMagnitude - 1;%限幅
            M_v2c_shifted(M_v2c_shifted < -maxMagnitude) = -maxMagnitude;

            % CNUs
            sorted_abs_values = sort(abs(M_v2c_shifted));
            min0 = sorted_abs_values(1, :);
            min1 = sorted_abs_values(2, :);
            sgn = prod(signStrict(M_v2c_shifted));
            [~, idx0] = min(abs(M_v2c_shifted)); % used for comparing

            M_c2v = zeros(size(M_v2c_shifted, 1), size(M_v2c_shifted, 2));
            M_c2v_shifted = zeros(size(M_v2c_shifted, 1), size(M_v2c_shifted, 2));
            for i = 1 : size(M_v2c_shifted, 1)
                for j = 1 : size(M_v2c_shifted, 2)
                    tmp = M_v2c_shifted(i, j);
                    if i ~= idx0(j) % instead of comparing two data
                        M_c2v(i, j) = min0(j) * sgn(j) * signStrict(tmp);
                    else
                        M_c2v(i, j) = min1(j) * sgn(j) * signStrict(tmp);
                    end
                end
                M_c2v_shifted(i, :) = circshift(M_c2v(i, :), [0, HBG( l, read_col_idx(i) )]);
            end

            % Write to llr_ram and M_c2v_l_it
            for i = 1 : length(read_col_idx)
                col = read_col_idx(i);
                start_idx = Zc * (col-1) + 1;
                end_idx = Zc * col;
                M_c2v_ram(l, start_idx:end_idx) = M_c2v_shifted(i, :);
                llr_ram(start_idx:end_idx) = M_v2c_fifo(i, :) + M_c2v_shifted(i, :);
            end

%             M_c2v_ram
%             llr_ram

        end % all layers end
        
        % Soft decision
        signbit = zeros(1,C_HBG * Zc);
        signbit(llr_ram >= 0) = 0;
        signbit(llr_ram <  0) = 1;

%         syndrome = signbit * H';
%         mean(syndrome)
%         if mean(syndrome) < 2
%             break;
%         end


    end % itera

    y = signbit(1 : kb*Zc);
    llr_final = llr_ram(1:kb*Zc);

end