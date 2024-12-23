function y = cyclic_shifter(llr, shift_ctrl, Zc, N)
    llr_shifted = zeros(1, length(llr));
    for i = 1 : N
        start_idx = (i-1) * Zc + 1;
        end_idx   = i * Zc;

        llr_group = llr(start_idx : end_idx);

        shift_amount = shift_ctrl(i);

        llr_shifted(start_idx : end_idx) = circshift(llr_group, [0, shift_amount]);
    end
    y = llr_shifted;
end