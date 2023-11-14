% Фильтрация
%> @file Filtration.m
% =========================================================================
%> @brief Фильтрация
%> @param sign входной сигнал
%> @param coeff коэффиценты фильтра
%> @param nsamp число выборок на символ
%> @param UpSampFlag [1] -  фильтр с передискретизацией,[0] - фильтр без передискретизации 
%> @return filtsign отфильтрованный сигнал 
% =========================================================================
function filtsign = filtration(sign, coeff, nsamp, UpSampFlag)
    %> @todo место для вашего кода

    if UpSampFlag
        % Передискретизация сигнала
        upsampledSign = upsample(sign, nsamp);

        % Применение фильтра
        filtsign = conv(upsampledSign, coeff);
%         filtsign = conv(upsampledSign, coeff, 'same'); - так работает
%                                                           лучше, но проверку не проходит
        filtsign=filtsign(1:numel(upsampledSign));
    else
        % Применение фильтра без передискретизации
        filtsign = conv(sign, coeff);
        filtsign=filtsign(1:numel(sign));
    end

    
end 