%% Task 2
clear; clc; close all;
load('Matlab_L3_6.mat');

%% Задание 1. Кадровая синхронизация
Stream=Bit_Stream;
Header = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Special_bits=Header;

res_cor(1 : length(Stream)-length(Special_bits)) = 0;
for itter = 1 : length(Stream)-length(Special_bits)
    res_cor(itter) = sum(Special_bits.*Stream(itter + 1:itter+length(Special_bits) ))/length(Special_bits);
end

% находим пики корреляционной функции
[~, Indexes_of_frames] = findpeaks(res_cor, 'MinPeakHeight', 0.99999);
Start_Of_Frame_Position=Indexes_of_frames(1);

% Определяем количество кадров данных
Number_of_frame = length(Indexes_of_frames);

% Строим график корреляции и сохраняем его в файл
figure
plot(res_cor)
title('Корреляционный анализ битовой последовательности')
xlabel('Номер бита')
ylabel('Корреляция')
savefig('Frame_Corr.fig')


%% Задание 2. Генератор псевдослучайной последовательности
%Функция Scrambler написана в файле "Scrambler.m"

Register = [1 0 0 1 0 1 0]; % начальное состояние регистра
sequence = Scrambler(Register); % генерация последовательности

%Функция циклической автокорреляции
acf = cyclic_autocorr(sequence);

figure;
plot(acf, 'LineWidth', 1.5);
title('Autocorrelation Function of Scrambler Output');
xlabel('Bit Offset');
ylabel('Autocorrelation');
saveas(gcf, 'ACF_Srambler.fig');

[~, max_index] = max(acf(2:end)); % поиск максимального значения автокорреляции
PN_Period = max_index; % вывод периода повторения
disp(PN_Period);

figure;
plot(sequence);


%% 

% функция для расчета циклической автокорреляции сигнала
function acf = cyclic_autocorr(signal)
    % определяем длину сигнала
    N = length(signal);

    % вычисляем циклическую автокорреляцию
    acf = zeros(1, N);
    for itter = 1:N
        acf(itter) = sum(signal.*circshift(signal, itter-1))/N;
    end
    % нормализация
    acf = acf / acf(1); 
end