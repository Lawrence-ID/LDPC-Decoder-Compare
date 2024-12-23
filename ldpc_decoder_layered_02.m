function [y, llr_final, it] = ldpc_decoder_layered_02(llr_data, llr_width, intera, HBG, H, check_ctrl, Zc)
    llr_ram = llr_data(:)'; % Init llr_ram

%     HBG=[1 1 1 1 -1 0 1; 0 0 -1 0 0 0 0;2 -1 2 2 2 2 2]; Zc = 3;kb=4;LDPC_K = kb * Zc;LDPC_Q = 3;LDPC_N = Zc * (kb + LDPC_Q);
%     [H, ~, ~] = ldpc_02_HMatrix_calc(HBG, Zc);LDPC_P = Zc * LDPC_Q;
%     llr_ram = [-10,-34,45,16,-34,30,-60,2,-10,20,-17,-39,12,-5,-14,-7,11,44,4,8,43];

%     HBG=[1 -1 2 1; 0 2 -1 0;2 -1 1 1]; Zc = 3; kb=1; LDPC_K = kb * Zc;LDPC_Q = 3;LDPC_N = Zc * (kb + LDPC_Q);
%     check_ctrl = [3, 3, 3];
%     [H, ~, ~] = ldpc_02_HMatrix_calc(HBG, Zc);
%     llr_ram = [-10,-34,45,16,-34,30,-60,2,-10,20,-17,-39];
%     intera = 2;

    % HBG: R*C, H: M*N = (R*Zc)*(C*Zc)
    [R_HBG, C_HBG] = size(HBG);
    kb = C_HBG - R_HBG;

    % Init M_c2v_ram: Size = 46 * Zc_max * 16[min0|min1|idx0|gsgn] bits
    M_c2v_ram = cell(R_HBG, Zc);
    for i = 1 : R_HBG
        for j = 1 : Zc
            M_c2v_ram{i, j} = zeros(1, 4);
        end
    end
    
    % Init M_c2v_sign_ram: Size = sum(check_ctrl) * Zc_max bits (sum(check_ctrl) is the total num of all edges)
    M_c2v_sign_ram = ones(R_HBG * Zc, C_HBG * Zc);

    % V2C RAM: Size = max(check_ctrl) * (8+1)bits * Zc_max 
    M_v2c_ram = zeros(max(check_ctrl), Zc);
    M_v2c_vaddr_ram = zeros(max(check_ctrl), Zc);

    maxMagnitude = 2^(llr_width-1);

    for it = 1 : intera
%         it
        for l = 1 : R_HBG

            tp = 1;

            new_min0 = zeros(1, Zc);
            new_min1 = zeros(1, Zc);
            new_idx0 = zeros(1, Zc);
            new_global_sign = zeros(1, Zc);

            for n = 1 : C_HBG % get Zc LLR at once to VNUs
                if HBG(l, n) ~= -1
                    start_idx = Zc * (n-1) + 1;
                    end_idx = Zc * n;
                    llr_shifted = circshift(llr_ram(start_idx:end_idx), [0, -HBG(l, n)]); % shift LLR first

                    cnt = HBG(l, n);

                    for j = 1 : Zc
                        i = start_idx + cnt;
                        old_min0 = M_c2v_ram{l, j}(1); old_min1 = M_c2v_ram{l, j}(2); old_idx0 = M_c2v_ram{l, j}(3); old_global_sign = M_c2v_ram{l, j}(4); % M_c2v_l_it-1, 用到上一次迭代的数据
                        M_c2v_sign_old = M_c2v_sign_ram((l-1) * Zc + j, i); % edge: cj<->vi
                        if i == old_idx0
                            M_ci_vj_old = old_global_sign * M_c2v_sign_old * old_min1;
                        else
                            M_ci_vj_old = old_global_sign * M_c2v_sign_old * old_min0;
                        end
                        cnt = mod(cnt+1, Zc);

                        % VNU, Calc M_vi_cj(9bit) = LLR_vi(8bit) - M_ci_vj_old(6bit)
                        M_vi_cj = llr_shifted(j) - M_ci_vj_old;
                        

                        % Store in V2C RAM
                        M_v2c_ram(tp, j) = M_vi_cj; % 第l行第tp个非-1元素对应着Zc个vi
                        M_v2c_vaddr_ram(tp, j) = i;
                       
                        abs_val = abs(M_vi_cj);     % 8bit
                        abs_val = max(min(abs_val, maxMagnitude-1), 0); % saturate
                        sign_val = signStrict(M_vi_cj);

                        % CNU, Update M_c2v
                        if tp == 1
                            new_min0(j) = abs_val;
                            new_min1(j) = abs_val;
                            new_idx0(j) = i;
                            new_global_sign(j) = sign_val;
                        else
                            if abs_val < new_min0(j)
                                new_min1(j) = new_min0(j);
                                new_min0(j) = abs_val;
                                new_idx0(j) = i;
                            elseif abs_val < new_min1(j) || tp == 2
                                new_min1(j) = abs_val;
                            end
                            new_global_sign(j) = new_global_sign(j) * sign_val;
                        end

%                         fprintf("M_c2v_sign_ram(c%d, v%d) = %d\n", (l-1)*Zc + j, i, M_vi_cj);
                        M_c2v_sign_ram((l-1) * Zc + j, i) = sign_val;
                        if tp == check_ctrl(l)
%                             new_min0_sat = max(min(new_min0(j), 15), -15); % saturate to 5bit
%                             new_min1_sat = max(min(new_min1(j), 15), -15); % saturate to 5bit
                            M_c2v_ram{l, j} = [new_min0(j), new_min1(j), new_idx0(j), new_global_sign(j)];
                            % 与c((l-1)*Zc + j)连接的Zc个vi: 更新LLR_vi
                            for dc_cnt = 1 : check_ctrl(l)
                                M_v2c_vaddr_delayed = M_v2c_vaddr_ram(dc_cnt, j);
                                M_v2c_delayed       = M_v2c_ram(dc_cnt, j);
                                if M_v2c_vaddr_delayed == new_idx0(j)
                                    M_cj_vi = new_min1(j) * new_global_sign(j) * signStrict(M_v2c_delayed);
                                else
                                    M_cj_vi = new_min0(j) * new_global_sign(j) * signStrict(M_v2c_delayed);
                                end
                                llr_ram(M_v2c_vaddr_delayed) = M_v2c_delayed + M_cj_vi; % reshifted here by M_v2c_vaddr_delayed
                            end
                        end

                    end % cj

                    tp = tp + 1;
                end % HBG(l, n) ~= -1
            end % n = 1 : C_HBG

%         M_c2v_ram        
%         llr_ram

        end % all layers end
        
        % Soft decision
        signbit = zeros(1,C_HBG * Zc);
        signbit(llr_ram >= 0) = 0;
        signbit(llr_ram <  0) = 1;


    end % itera

    y = signbit(1 : kb*Zc);
    llr_final = llr_ram(1:kb*Zc);

end