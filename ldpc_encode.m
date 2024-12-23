% 5g ldpc encoding
% input: s, bit sequence of dimension K * 1
% output: encoded_bits, bit sequence
% reference: 3GPP TS 38.212 section 5.3.2
% author: Xiao, Shaoning ������
% license: MIT

%�������룺����������������s�������ͼ����
%�������������������s��У�����H��ZC��ֵ�������Ķ���������(s,w)

function [encoded_bits, H, Z_c, encoded_bits_original] = ldpc_encode(s, base_graph_index,kb)

K = length(s);  %������������еĳ���

encoded_bits = zeros(3*K, 1);%����һ��������������Ϊ3K
%BG1�Ĵ�СΪ46*68

if base_graph_index == 1    %�����ͼѡ��1
    a = 4;                  %BG1H�ĺ��ľ��󲿷�Ϊ4*22
    b = 22;                 %BG1�Ŀɱ���ı���Ϊ22Zc
    c = 26;
    d = 42;
    e = 46;
    z = K/b;                %Zc����ֵ����׼�ļ��㷽���������ֵ��Ҫ�����������Zc��ֵ���ܲ��Ǳ�׼��������ģ���Ҫ��������
    Z_c = z;
    N = 66 * Z_c;
    z = K/b;
    set_index = lifting_size_table_lookup(z);%��Ҫ������ILS�������ڱ�׼�ж����˰���ILS����������ʵ�ʵ�BG1���а�����ʽ�����������ʵ��ʵ�ַ����Ƚ����棬����ȱ���ǵ���֧�ֱ�׼��������ļ���Zֵ��ILS��������������ֵ�����ϱ�׼����ֵ��ô����
    load parity_check_matrices_protocol_1    %�ļ��д洢��BG1���л�ͼ��ʽ����������ȷ��ֵ��ΧΪ1~8
    BG = parity_check_matrices_protocol_1(:, :, set_index); %#ok<NODEF>   ѡ��BG1�Ļ�ͼ
%BG1�Ĵ�СΪ42*52

elseif base_graph_index == 2%          
    a = 4;                  %BG2H�ĺ��ľ��󲿷�Ϊ4*22
    b = 10;                 %BG2�Ŀɱ���ı���Ϊ10Zc
    c = 14;
    d = 38;
    e = 42;
    z = K/b;                %Zc����ֵΪ
    Z_c = z;
    N = 50 * Z_c;
    set_index = lifting_size_table_lookup(z); %��Ҫ������ILS�������ڱ�׼�ж����˰���ILS����������ʵ�ʵ�BG2���а�����ʽ   Э��38 212 19ҳ
    load parity_check_matrices_protocol_2     %�ļ��д洢��BG2���л�ͼ��ʽ����������ȷ��ֵ��ΧΪ1~8
    BG = parity_check_matrices_protocol_2(:, :, set_index); %#ok<NODEF>     ѡ��BG2�Ļ�ͼ
else
  error('wrong base graph index in ldpc encoding.');
end

BG(BG ~= -1) = mod(BG(BG ~= -1), Z_c); %ɧ�������ǲ��ǿ���������ʵ�ܼ򵥣����ǰ�ѡ�õĻ�ͼ��Zcȡ�࣬�����BG����Ӧ��ѭ����λ����
%�������ܴ��BGǰ���У����Ұ�-1Ԫ��ת��Ϊ0Ԫ��---��˼��BG��ǰ���д��������---����ʵ�ʱ��뻹���õ�


% set_index = lifting_size_table_lookup(Z_c);
% BG = parity_check_matrices_protocol(:, :, set_index);
%�����ǰ�BG������зֿ鴦��
%����Ļ�����ʽ�����ļ���ͼƬ
A_prime = BG(1:a, 1:b);             %����A--���ľ��󲿷�    a*b
B_prime = BG(1:a, (b+1):c);         %����B--˫�Խǽṹ����      a*a
C_prime = BG((a+1):e, 1:b);         %����C--��չ����a single parity-check��
D_prime = BG((a+1):e, (b+1):c);     %����D--��չ����a single parity-check��  C+D

if base_graph_index == 2 %��bg2ʱ�޸�bֵ
    if(kb == 10)        
        b = 10;
        s = s;
        K = K;
    elseif(kb == 9)
        b = 9;
        s = s(1:9*Z_c);
        K = K - 1 * Z_c;
    elseif(kb == 8)
        b = 8;
        s = s(1:8*Z_c);
        K = K - 2 * Z_c;
    elseif(kb == 6)
        b = 6;
        s = s(1:6*Z_c);
        K = K - 4 * Z_c;
    else
        error('wrong kb in ldpc encoding.');
    end
end

%���������kbȡֵ����������A��C��ȡֵ��Χ
A_prime = A_prime(:, 1:b);             %����A--���ľ��󲿷�    a*b
C_prime = C_prime(:, 1:b);         %����C--��չ����a single parity-check��

encoded_bits(1:K-2*Z_c) = s(2*Z_c+1:end);
% for k = (2*Z_c):(K-1)
%   if s(k+1) ~= -1
%     encoded_bits(k-2*Z_c+1) = s(k+1);
%   else
%     s(k+1) = 0;
%     encoded_bits(k-2*Z_c+1) = -1;
%   end    
% end

z = Z_c;
%spalloc�������ܣ�Ϊϡ��������ռ�--spalloc(10,10,20)---����Ԫ�����Ϊ20
%A�Ĵ�СΪ--a*z��b*z-
%size(A_prime)����a*z��b*z
%ones(size(A_prime))---����һ��a*z��b*z��С�ĵ�λ����
%Ϊʲô��ô�ҳ���������з������ֵ??
A = spalloc(a*z, b*z, nnz(A_prime + ones(size(A_prime))));
for row_index = 1:a     %������  1:a
    for column_index = 1:b      %������     1:b
        if A_prime(row_index, column_index) ~= -1  %if(�����ڲ�����ֵ������-1)  begin---����Ĺ����ǽ�A�����еķ���Ԫ��--���ݱ�׼�����---ͨ��ѭ����λ---�õ���չ��ľ���
            A((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(A_prime(row_index, column_index), z)+1):z, 1:mod(A_prime(row_index, column_index), z)], ones(1, z), z, z);
        end
    end
end
%����B�������4*4��С
B = spalloc(a*z, a*z, nnz(B_prime + ones(size(B_prime))));

for row_index = 1:a
    for column_index = 1:a
        if B_prime(row_index, column_index) ~= -1
            B((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(B_prime(row_index, column_index), z)+1):z, 1:mod(B_prime(row_index, column_index), z)], ones(1, z), z, z);
        end
    end
end
% d=e-a  ʣ�µ�����
C = spalloc(d*z, b*z, nnz(C_prime + ones(size(C_prime))));

for row_index = 1:d
    for column_index = 1:b
        if C_prime(row_index, column_index) ~= -1
            C((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(C_prime(row_index, column_index), z)+1):z, 1:mod(C_prime(row_index, column_index), z)], ones(1, z), z, z);
        end
    end
end

D = spalloc(d*z, a*z, nnz(D_prime + ones(size(D_prime))));

for row_index = 1:d
    for column_index = 1:a
        if D_prime(row_index, column_index) ~= -1
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = sparse(1:z, [(mod(D_prime(row_index, column_index), z)+1):z, 1:mod(D_prime(row_index, column_index), z)], ones(1, z), z, z);
        else
            D((row_index-1)*z+1:row_index*z, (column_index-1)*z+1:column_index*z) = spalloc(z, z, 0);
        end
    end
end

B_inv = spalloc(a*z, a*z, 20*z);

%speye�������ܣ�����һ��az*az�ĶԽǾ�������Ǿ������任B����任
%����һ�����  �������ļ�
%ÿ����ͼ��Ӧ��4*4����˫�Խǽṹ�ľ����Ԫ�ؾ�����ͬ���ⲿ����Ҫ���⴦��
if (base_graph_index == 1) && (set_index ~= 7)    %���ѡ��BG1---set_index��=7---BG1��ͼ��˫�Խǽṹ���������֣�����ֿ�����
    B_inv(1:z, 1:z)             = speye(z);       %B_inv�����z*zλ��Ϊ��λ����
    B_inv(1:z, 1+z:2*z)         = speye(z);       %����z*z�ĶԽǾ������  ��СΪz*2z
    B_inv(1:z, 1+2*z:3*z)       = speye(z);       %����z*z�ĶԽǾ������  ��СΪz*3z
    B_inv(1:z, 1+3*z:4*z)       = speye(z);       %�ĸ�z*z�ĶԽǾ������  ��СΪz*4z
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
%��λ��Ļ����������ƶ�z-1λ
elseif (base_graph_index == 2) && ((set_index ~= 4) && (set_index ~= 8))
    B_inv(1:z, 1:z)             = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+z:2*z)         = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+2*z:3*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+3*z:4*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
elseif (base_graph_index == 1) && (z == 208)
    B_inv(1:z, 1:z)             = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+z:2*z)         = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+2*z:3*z)       = sparse(circshift(eye(208), 105));
    B_inv(1:z, 1+3*z:4*z)       = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1:z)         = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1:z)       = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(208) + sparse(circshift(eye(208), 105));
    B_inv(1+3*z:4*z, 1:z)       = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(circshift(eye(208), 105)); 
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(208) + sparse(circshift(eye(208), 105));
elseif (base_graph_index == 1) && ((z ~= 208) && (set_index == 7))
    B_inv(1:z, 1:z)             = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+z:2*z)         = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+2*z:3*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1:z, 1+3*z:4*z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [z, 1:(z-1)], ones(1, z), z, z);
elseif (base_graph_index == 2) && ((set_index == 4) || (set_index == 8))    
    B_inv(1:z, 1:z)             = speye(z);
    B_inv(1:z, 1+z:2*z)         = speye(z);
    B_inv(1:z, 1+2*z:3*z)       = speye(z);
    B_inv(1:z, 1+3*z:4*z)       = speye(z);
    B_inv(1+z:2*z, 1:z)         = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+z:2*z)     = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+2*z:3*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+z:2*z, 1+3*z:4*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1:z)       = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+z:2*z)   = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+2*z:3*z, 1+3*z:4*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1:z)       = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+z:2*z)   = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+2*z:3*z) = sparse(1:z, [2:z, 1], ones(1, z), z, z);
    B_inv(1+3*z:4*z, 1+3*z:4*z) = speye(z) + sparse(1:z, [2:z, 1], ones(1, z), z, z);    
end

s = s(:);
%��ʽ---���������Ǿ������ķ���   B_invΪB�������---
p_1 = mod(B_inv * (A * s), 2);
p_2 = mod(C * s + D * p_1, 2);
w = [p_1; p_2];

%�ѱ���ı��ز�ȫУ��λ
% for k = K:(N+2*Z_c-1)
%    encoded_bits(k-2*Z_c+1) = w(k-K+1);   
% end    
encoded_bits(K-2*Z_c+1:K-2*Z_c+length(w)) = w;
%����H����  ǰ���Ѿ���ȥ��ǰ���е�BG---BG1 46*66  BG2---42*52�����ɵ�H����
%spalloc(a*z, d*z, 0)����0��        speye(d*z)���ɶԽǾ���
H = [A, B, spalloc(a*z, d*z, 0); C, D, speye(d*z)];

%ԭʼ�ı�����أ�BG��Ӧ��ǰ���б��������
encoded_bits_original = [s; w]; 

% mod(H * [s; p_1; p_2]) = 0

clear A
clear B
clear C
clear D
clear B_inv
clear A_prime
clear B_prime
clear C_prime
clear D_prime

end

%ǰ��������Ĳ���  Z=4�����ɽ������ͼ��ʾ��
%1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     0     1     0
%0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     1     1     0     0     0     1     0     0     0     1     0     0     0     1     0
%0     0     1     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     0     0     1     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     0     0     0     1     0     0     1     1     0     0     1     1     0     0
%0     0     1     0     0     0     1     0     0     1     1     0     0     1     1     0
%0     0     0     1     0     0     0     1     0     0     1     1     0     0     1     1
%1     0     0     0     1     0     0     0     1     0     0     1     1     0     0     1
%0     1     0     0     0     1     0     0     0     1     0     0     1     1     0     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     1     1     0
%0     0     0     1     0     0     0     1     0     0     0     1     0     0     1     1
%1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     1


%�����ƶ�3λ
%0     0     0     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     0     1     0
%1     0     0     1     0     0     0     1     0     0     0     1     0     0     0     1
%1     1     0     0     1     0     0     0     1     0     0     0     1     0     0     0
%0     1     1     0     0     1     0     0     0     1     0     0     0     1     0     0
%0     0     1     1     0     0     1     0     0     0     1     0     0     0     1     0
%1     0     0     1     1     0     0     1     0     0     0     1     0     0     0     1
%1     1     0     0     1     1     0     0     1     0     0     0     1     0     0     0
%0     1     1     0     0     1     1     0     0     1     0     0     0     1     0     0
%0     0     1     1     0     0     1     1     0     0     1     0     0     0     1     0
%0     0     0     1     0     0     0     1     0     0     0     1     1     0     0     1
%1     0     0     0     1     0     0     0     1     0     0     0     1     1     0     0
%0     1     0     0     0     1     0     0     0     1     0     0     0     1     1     0
%0     0     1     0     0     0     1     0     0     0     1     0     0     0     1     1
