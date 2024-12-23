% 5g ldpc parameter calculation
% input: Z
% output: set_index
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning 萧少宁
% license: MIT





%函数功能：计算出ILS索引
%输入参数：根据码块长度计算出的Z值
function set_index = lifting_size_table_lookup(Z)

first_element = bitshift(Z, -ntz(Z));%bitshift函数功能：如bitshift(2,1)=4---2的平方

lut = [1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8];

if first_element > 15
  error('wrong input parameter Z in lift_size_table_lookup.');    
end

set_index = lut(first_element);

if set_index == 0
  error('wrong input parameter Z in lift_size_table_lookup.');    
end

end