% Генерация коэффицентов (импульсной характеристики)для фильтра корень 
% из приподнятого косинуса
%> @file sqRCcoeff.m
% =========================================================================
%> @brief Генерация коэффицентов (импульсной характеристики)для фильтра корень 
%> из приподнятого косинуса
%> @param span Длина фильтра в символах (количество боковых лепестков sinc, 
%> сумма с двух сторон)
%> @param nsamp Число выборок на символ
%> @param rolloff Коэффицент сглаживания (alfa)
%> @return coeff коэффиценты для фильтра корень из приподнятого косинуса
% =========================================================================
function coeff = sqRCcoeff (span, nsamp, rolloff)
    Ts=1;
    % Вектор времени    
    t1 = (-span/2):(1/nsamp):(span/2);
    coeff = zeros(1, numel(t1));

    % Расчет коэффициентов корня из приподнятого косинуса
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

    % Нормализация коэффициентов фильтра
    coeff = coeff / sqrt(sum(coeff.^2));

end