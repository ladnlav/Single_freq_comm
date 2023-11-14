% ����������
%> @file Filtration.m
% =========================================================================
%> @brief ����������
%> @param sign ������� ������
%> @param coeff ����������� �������
%> @param nsamp ����� ������� �� ������
%> @param UpSampFlag [1] -  ������ � ������������������,[0] - ������ ��� ����������������� 
%> @return filtsign ��������������� ������ 
% =========================================================================
function filtsign = filtration(sign, coeff, nsamp, UpSampFlag)
    %> @todo ����� ��� ������ ����

    if UpSampFlag
        % ����������������� �������
        upsampledSign = upsample(sign, nsamp);

        % ���������� �������
        filtsign = conv(upsampledSign, coeff);
%         filtsign = conv(upsampledSign, coeff, 'same'); - ��� ��������
%                                                           �����, �� �������� �� ��������
        filtsign=filtsign(1:numel(upsampledSign));
    else
        % ���������� ������� ��� �����������������
        filtsign = conv(sign, coeff);
        filtsign=filtsign(1:numel(sign));
    end

    
end 