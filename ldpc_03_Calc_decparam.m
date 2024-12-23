function [interweave_addr,barrel_shifter,variable_ctrl,check_ctrl, V] = ldpc_03_Calc_decparam(LDPC_N,LDPC_K,LDPC_P,LDPC_Q,H_matrix, Zc)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Parameters
    % LDPC_N                    :           Code Length
    % LDPC_K                    :           Information Bits Length
    % LDPC_P                    :           Parity Check Bis
    % LDPC_Q                    :           Rows of BG
    % H_matrix                  :           Parity Check Matrix
    % Zc                        :           Expension Factors
    % Output Parameters     
    % check_ctrl                :           Number of columns to store log domain messages.
    % variable_ctrl             :           buf_d = [40	124	281	374	717	950	1041 1145 19 340 420 623 706 746 847 1013 1145 1249]
    %                           :           Zc = 104
    %                           :           variable_ctrl = [2 1 ...]
    %                           :           variable_ctrl is the sum of non-zero elements in each group located in different Zc regions.
    % interweave_addr           :           interweave_addr is the index of the non-zero element in each group located in a different Zc region.
    % shuffer_ctrl_s            :           Shift coefficient of barrel shift register.
    %                           :           For variable node updates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate check_ctrl
    check_ctrl = zeros(1,LDPC_Q);
    for i = 1:LDPC_Q
        check_ctrl(i)   =   sum(H_matrix((i-1)*Zc+1,:));
    end
    % group_len: The number of log field messages per group.
    group_len           =   sum(check_ctrl);
    %% Break and reorganize the parity check matrix and divide it into Zc groups.
    H_T                 =   zeros(LDPC_P,LDPC_N);
    % For example: H_T(1) = H_matrix(1)    H_T(2) = H_matrix(1+Zc)     ......
    for mm = 1 : Zc
        for nn = 1 : LDPC_Q
            H_T((mm-1)*LDPC_Q+nn,:) = H_matrix((nn-1)*Zc+mm,:);
        end
    end
    %% Find the index of the non-zero element in H_T, stored in buf_a.
    buf_a               =   zeros(LDPC_P,LDPC_N);
    for i = 1:LDPC_P
        index_a = find(H_T(i,:)==1);
        buf_a(i,index_a) = index_a;
    end
    %% Merge each group in buf_a into one row.
    buf_b = (reshape(buf_a',LDPC_N*LDPC_Q,Zc))';
    %% Calculate the index of the non-zero element in buf_b.
    buf_c = zeros(Zc,group_len);
    for j = 1:Zc
        index_b = find(buf_b(j,:) ~= 0);
        buf_c(j,:) = buf_b(j,index_b);
    end
    %% Get a row in buf_c
    buf_d = buf_c(1,:);
    %% control word calculation
    variable_ctrl_s     = zeros(1,LDPC_N/Zc);
    inter_wv_addr_s     = zeros(1,length(buf_d));
    shuffer_ctrl_s      = zeros(1,length(buf_d));
    base = 0;
    for ii = 1:LDPC_N/Zc
        index = find((buf_d <= Zc*ii) & (buf_d >= Zc*(ii-1)+1));         
        variable_ctrl_s(ii) = length(index);
        for jj = 1:variable_ctrl_s(ii)
            inter_wv_addr_s(base+jj) = index(jj)-1;
            shuffer_ctrl_s(base+jj) = mod(buf_d(index(jj)),Zc)-1;
            temp_shuffer = mod(buf_d(index(jj)),Zc);
            if(temp_shuffer>=0)
                shuffer_ctrl_s(base+jj) = temp_shuffer;
            else
                shuffer_ctrl_s(base+jj) = Zc+temp_shuffer;
            end
        end
        base = base + variable_ctrl_s(ii);
    end
    variable_ctrl = variable_ctrl_s ;
    interweave_addr = inter_wv_addr_s;
    barrel_shifter = shuffer_ctrl_s;
    %% Get V(ci): Set of variable nodes connected with ci
    V = cell(size(H_matrix, 1), 1);
    for i = 1 : size(H_matrix, 1)
        indices = find(H_matrix(i, :) == 1);
        V{i} = indices;
    end
end
