%% Frame synchronization lab
clear; clc; close all;

%% Task 1: Autocorrelation of given m-seq

Register = [0 1 0 0 1 1 0]; % начальное состояние регистра
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

% figure;
% plot(sequence);

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