function sequence = Scrambler(Register)
    % Инициализация параметров
    %m = 8; % число бит регистра
    sequence_length = 2^7 - 1; % длина псевдослучайной последовательности
    sequence = zeros(1, sequence_length); % инициализация последовательности
    
    % Генерация псевдослучайной последовательности
    for i = 1:sequence_length
        sequence(i)=array_xor(Register); % выходной бит
        feedback = sequence(i); 

        Register = circshift(Register,1); % сдвиг регистра
        Register(1) = feedback; % обновление первого бита
    end
end

function xor_result = array_xor(Register)
    
    coeff_mask=[1 1 1 0 1 1 1 1]; %367 (oct) in bits
    xor_in=Register(logical(coeff_mask(2:end)));

    xor_result = xor_in(1);
    for i = 2:length(xor_in)
        xor_result=xor(xor_result,xor_in(i));
    end

end