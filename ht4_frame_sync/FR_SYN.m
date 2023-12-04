%% Frame synchronization lab
clear; clc; close all;

%% Task 1: Autocorrelation of given m-seq

Register = [0 1 0 0 1 1 0]; % начальное состояние регистра
sequence = Scrambler(Register); % генерация последовательности

%Функция циклической автокорреляции
acf = cyclic_autocorr(sequence);

% figure;
% plot(acf, 'LineWidth', 1.5);
% title('Autocorrelation Function of Scrambler Output');
% xlabel('Bit Offset');
% ylabel('Autocorrelation');
% saveas(gcf, 'ACF_Srambler.fig');

[~, max_index] = max(acf(2:end)); % поиск максимального значения автокорреляции
PN_Period = max_index; % вывод периода повторения
disp(PN_Period);

% figure;
% plot(sequence);

%% Task 2: Signal generation + Probability of detection
Header=sequence;
Amount_of_Frame = 1001; % one frame for drop
L_H=length(Header);    % length of header
Length_Data_IQ = 10*length(Header); % length of data
L_D=L_H+Length_Data_IQ;         %frame length



Tx_Bits = randi([0 1], 1, Amount_of_Frame*Length_Data_IQ); % генерация бит; 

TX_IQ_Data = mapping(Tx_Bits, 'BPSK');
Header_IQ = mapping (Header,'BPSK');
% Frame structure 
% |L_H Header| L_D=10*L_H Data|
IQ_TX_Frame = FrameStruct(TX_IQ_Data, Header_IQ, Amount_of_Frame);
IQ_TX_Frame = IQ_TX_Frame(L_D-30:end); % dropping first frame, but saving
                                        %  first 30 symbols for some offset

True_Indexes=31:L_D:L_D*(Amount_of_Frame-1); % setting true starts of the frames

SNR = -10:1:25;
prob_det=zeros(size(SNR));
Freq_Offset=0.1;
Rx_amount_frames=floor(length(IQ_TX_Frame)/L_D);

parfor i=1:length(SNR)
    snr=SNR(i);
    Indexes_of_frames = zeros([1,Rx_amount_frames]);

    % Channel
    Channel_IQ = awgn(IQ_TX_Frame, snr, 'measured');
    %Channel_IQ = Channel_IQ.*exp(1j*2*(1:size(Channel_IQ,2))*pi*Freq_Offset);
    
    cross_corr=corr(Header_IQ,Channel_IQ);
    
    % finding possible starts of the frames
    for k=1:Rx_amount_frames

        if k==Rx_amount_frames
            [~, index_i] = max(cross_corr(1+(k-1)*L_D:end));
            Indexes_of_frames(k)=(k-1)*L_D+index_i;
            continue;
        end

        [~, index_i] = max(cross_corr(1+(k-1)*L_D:k*L_D));
        Indexes_of_frames(k)=(k-1)*L_D+index_i;
    end

    % calculating probability of right detection
    num_true=sum(Indexes_of_frames==True_Indexes);
    prob_det(i)=num_true/length(True_Indexes);

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



