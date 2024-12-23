%% 5G NR LDPC���� (֧��ȫ����)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������
%   srcData  ����Դ����
%   BG1or2   ��ѡ���ͼ��1��BG1,2��BG2
%   Z        ���������ӣ����жȣ�
% �������
%   HBG      ��У�������
%   K        ����Ϣλ����

%HBG: base matrix
%Z: expansion factor
%srcData: message vector, length = (#cols(HBG)-#rows(HBG))*Z
%cword: codeword vector, length = #cols(HBG)*Z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cword = ldpc_04_Encoder(HBG,Z,srcData)
    [m,n] = size(HBG);

    cword = zeros(1,n*Z);
    cword(1:(n-m)*Z) = srcData;

    %double-diagonal encoding
    % calc AS
    temp = zeros(1,Z);
    for i = 1:4 %row 1 to 4
        for j = 1:n-m %message columns
            temp = mod(temp + mul_sh(srcData((j-1)*Z+1:j*Z),HBG(i,j)),2);
        end
    end
    
    % B * As
    if HBG(2,n-m+1) == -1
        p1_sh = HBG(3,n-m+1);
    else
        p1_sh = HBG(2,n-m+1);
    end
    cword((n-m)*Z+1:(n-m+1)*Z) = mul_sh(temp,Z-p1_sh); %p1
    
    %Find p2, p3, p4
    for i = 1:3
        temp = zeros(1,Z);
        for j = 1:n-m+i
            temp = mod(temp + mul_sh(cword((j-1)*Z+1:j*Z),HBG(i,j)),2);
        end
        cword((n-m+i)*Z+1:(n-m+i+1)*Z) = temp;
    end
    
    %Remaining parities
    for i = 5:m
        temp = zeros(1,Z);
        for j = 1:n-m+4
            temp = mod(temp + mul_sh(cword((j-1)*Z+1:j*Z),HBG(i,j)),2);
        end
        cword((n-m+i-1)*Z+1:(n-m+i)*Z) = temp;
    end
end
%% ֱ����λ����
% x: input block
% k: -1 or shift
% y: output
function y = mul_sh(x,k)


    if (k==-1)
        y = zeros(1,length(x));
    else
        y = [x(k+1:end) x(1:k)]; %multiplication by shifted identity
    end
end
