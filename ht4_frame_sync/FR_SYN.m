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

%% Task 2: Signal generation + Probability of detection
Header=sequence;
Amount_of_Frame = 1001;
L_H=length(Header);
Length_Data_IQ = 10*length(Header);
L_D=L_H+Length_Data_IQ;

True_Indexes=L_D:L_D:L_D*(Amount_of_Frame-1);

Tx_Bits = randi([0 1], 1, Amount_of_Frame*Length_Data_IQ); % генерация бит; 

% Frame structure 
% |L_H Header| L_D=10*L_H Data|
IQ_TX_Frame = FrameStruct(Tx_Bits, Header, Amount_of_Frame);

SNR = -10:0.1:5;
prob_det=zeros(size(SNR));
Freq_Offset=0.1;

parfor i=1:length(SNR)
    snr=SNR(i);
    
    prob_det_m=zeros(size(1:10));
    for j=1:10
        % Channel
        Channel_IQ = awgn(IQ_TX_Frame, snr, 'measured');
        %Channel_IQ = Channel_IQ.*exp(1j*2*(1:size(Channel_IQ,2))*pi*Freq_Offset);
        
        cross_corr=corr(Header,Channel_IQ);
        
        % находим пики корреляционной функции
        [~, Indexes_of_frames] = findpeaks(cross_corr, 'SortStr', 'descend', 'NPeaks', length(True_Indexes  ));
        Indexes_of_frames = sort(Indexes_of_frames);
        
        num_true=sum(Indexes_of_frames==True_Indexes);
    
        prob_det_m(j)=num_true/length(True_Indexes);

    end
    prob_det(i)=mean(prob_det_m);
end

figure(3);
plot(SNR,prob_det, '-o','LineWidth', 2,'MarkerSize', 2);
xlabel('SNR (dB)');
ylabel('Probability of detection');
title('Probability of detection vs SNR');
grid on;
hold on;


%% Functions

function res_cor = corr(header, signal)
    res_cor(1 : length(signal)-length(header)) = 0;
    for itter = 1 : length(signal)-length(header)
        res_cor(itter) = sum(header.*signal(itter + 1:itter+length(header) ))/length(header);
    end

    res_cor=res_cor/max(res_cor);
end


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



