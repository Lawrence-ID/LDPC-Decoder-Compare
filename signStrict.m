function y = signStrict(x)
    y = sign(x);
    y(y == 0) = 1;  % 将0转换为1
end