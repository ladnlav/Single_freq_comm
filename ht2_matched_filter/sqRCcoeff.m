% ��������� ������������ (���������� ��������������)��� ������� ������ 
% �� ������������ ��������
%> @file sqRCcoeff.m
% =========================================================================
%> @brief ��������� ������������ (���������� ��������������)��� ������� ������ 
%> �� ������������ ��������
%> @param span ����� ������� � �������� (���������� ������� ��������� sinc, 
%> ����� � ���� ������)
%> @param nsamp ����� ������� �� ������
%> @param rolloff ���������� ����������� (alfa)
%> @return coeff ����������� ��� ������� ������ �� ������������ ��������
% =========================================================================
function coeff = sqRCcoeff (span, nsamp, rolloff)
    Ts=1;
    % ������ �������    
    t1 = (-span/2):(1/nsamp):(span/2);
    coeff = zeros(1, numel(t1));

    % ������ ������������� ����� �� ������������ ��������
    for idx = 1:length(t1)
        time = t1(idx);
        
        switch time
            case 0
                coeff(idx) = 1/Ts*(1+rolloff*(4/pi-1));
                
            case {Ts/(4*rolloff), -Ts/(4*rolloff)}
                coeff(idx)=rolloff/(Ts*sqrt(2))*((1+2/pi)*sin(pi/4/rolloff) ...
                    +(1-2/pi)*cos(pi/4/rolloff));
                
            otherwise
                coeff(idx) = 1/Ts * (sin(pi * time / Ts * (1-rolloff)) + ...
                    4 * rolloff * time / Ts * cos(pi * time / Ts * (1+rolloff))) / ...
                    (pi * time / Ts * (1 - (4 * rolloff * time / Ts)^2));
        end
    end

    % ������������ ������������� �������
    coeff = coeff / sqrt(sum(coeff.^2));

end