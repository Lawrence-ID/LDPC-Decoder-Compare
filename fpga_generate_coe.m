%% ���ݸ�����LDPC�����������ּ���FPGA��������ROM��coe�ļ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	�� �� �ˣ���  ��
%	�������ڣ�2021/05/26
%	�����־��
%       2021/05/26     �½�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����coe�ļ�
function fpga_generate_coe(rom_data,coe_name_str)
    file_name = strcat(coe_name_str,'.coe');
    fp = fopen(file_name,'wt');
    fprintf(fp,'memory_initialization_radix = 10;\nmemory_initialization_vector =\n');
    [M,N] = size(rom_data);
    for ii = 1 : M
        for jj = 1 : N
            fprintf(fp,'%d', rom_data(ii,jj));
        end
        if ii == M
            fprintf(fp,';\n');
        else
            fprintf(fp,',\n');
        end
    end
    fclose(fp);
end