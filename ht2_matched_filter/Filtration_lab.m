% ����������.
% =========================================================================
%> ���������� �������� �����
% =========================================================================
    %> �������� workspace
    clear all;
    %> �������� ��������
    close all;
    %> �������� Command Window
    clc;
% =========================================================================
%> ������� sinc. ������ �� ������.
% =========================================================================
    %> ������� ������
    x = [-10:0.1:10];
    %> ������� sinc
    y = sinc(x);
    % =====================================================================
    %> ������ ���������� ��������������
    % =====================================================================
    figure 
    plot(x,y)
    title('sinc(x)')
    % =====================================================================
    %> ������ �������
    % =====================================================================
    %> ����� + ������ � ��
    spectum = 10*log10(abs(fft(y)));
    %> ������ �� ������
    spectum = [spectum(102:201), spectum(1:101)];
    %> ������ � ������ T = 2 (��. �������� sinc)
    figure 
    plot(x/2, spectum)
    title('spectum sinc(x)')
% =========================================================================
%> ������ 1: �������� �������, ������� ���������� ����������� (���������� 
%> ��������������)��� ������� ������ �� ������������ ��������.
%> ��������� ���������� � ��������� �������������� �������.
% =========================================================================
    %> ����� ������� � �������� (����� ������� ��������� sinc, ����� � ���� ������)
    span = 20;
    %> ����� ������� �� ������
    nsamp  = 4;
    %> ���������� ����������� (alfa)
    rolloff = 0.2;
    % =====================================================================
    %> @todo ��������� �������
    sqimpuls = sqRCcoeff (span, nsamp, rolloff);
    %> @todo ��������� ���������� � ��������� �������������� �������

    % ���������� ��������������
    pulse = (sqimpuls);
    % ��������� �������������� �������
    frequency_response = abs(fftshift(fft(pulse)));
    
    figure;
    subplot(2,1,1);
    plot(pulse);
    title('Root Raised Cosine Pulse Response');
    xlabel('Sample Index');
    ylabel('Amplitude');
    
    subplot(2,1,2);
    f = linspace(-5, 5, length(frequency_response));
    plot(f, 10*log10(frequency_response));
    title('Frequency Response');
    xlabel('Frequency');
    ylabel('Magnitude');
    xlim([-5, 5]);
    ylim([-60, 10]);
% =========================================================================
%> �������� 1.
%> ��������� �� ����������� ��������
% =========================================================================
txfilter1 = comm.RaisedCosineTransmitFilter('RolloffFactor', rolloff, ...
                                           'FilterSpanInSymbols',span,...
                                           'OutputSamplesPerSymbol', nsamp);
check1 = coeffs(txfilter1);
if sum(abs(check1.Numerator-sqimpuls))< 0.001 % �������� ���������� ���� 
                                              % ���������� ������������� 
                                              % � �������� ���������
    ans = '�������� ������ 1 �������� �������'
else 
    err = '������ � ������ 1. ��������� �����������'
    ans = sum(abs(check1.Numerator-sqimpuls))
end
% =========================================================================
%> ������ 2: �������� �������, ������� ���������� ����������� (���������� 
%> ��������������) ��� ������� ������������ ��������.
%> ��������� ���������� � ��������� �������������� �������.
    %> ��������� ���������� �������������� ��� �����, ��� ����� � �������������� sinc 
    %> �� ����� �������
% =========================================================================
    %> ����� ������� � �������� (����� ������� ��������� sinc, ����� � ���� ������)
    span = 20;
    %> ����� ������� �� ������
    nsamp  = 4;
    %> ���������� ����������� (alfa)
    rolloff = 0.2;
    % =====================================================================
    %> @todo ��������� �������
    impuls = RCcoeff (span, nsamp, rolloff);
    %> @todo ��������� ���������� � ��������� �������������� �������

    % ���������� ��������������
    pulse2 = (impuls);
    % ��������� �������������� �������
    frequency_response2 = abs(fftshift(fft(pulse2)));
    
    figure;
    subplot(2,1,1);
    plot(pulse2);
    title('Raised Cosine Pulse Response');
    xlabel('Sample Index');
    ylabel('Amplitude');
    
    subplot(2,1,2);
    f = linspace(-5, 5, length(frequency_response2));
    plot(f, 10*log10(frequency_response2));
    title('Frequency Response');
    xlabel('Frequency');
    ylabel('Magnitude');
    xlim([-5, 5]);
    ylim([-60, 10]);

    t1 = (-span/2):(1/nsamp):(span/2);
    figure;
    plot(pulse);
    hold on;
    plot(pulse2);
    plot(sinc(t1*(1-rolloff)));
    title('All pulses in one chart');
    xlabel('Sample Index');
    ylabel('Amplitude');
    legend("Root Raised Cosine", "Raised Cosine", "Sinc",'Location','southwest');
% =========================================================================
%> �������� 2.
%> ��������� �� ����������� ��������
% =========================================================================
txfilter2 = comm.RaisedCosineTransmitFilter('RolloffFactor', rolloff, ...
                                            'FilterSpanInSymbols',span,...
                                            'OutputSamplesPerSymbol', nsamp,...
                                            'Shape', 'Normal');
check2 = coeffs(txfilter2);
if sum(abs(check2.Numerator-impuls))< 0.1 % �������� ���������� ���� 
                                          % ���������� ������������� 
                                          % � �������� ���������
    ans = '�������� ������ 2 �������� �������'
else 
    err = '������ � ������ 2. ��������� �����������'
    ans = sum(abs(check2.Numerator-impuls))
end
% =========================================================================
%> ������� 3. 
%> �������� ������� ����������, ������� �������� � ���� �������: �
%> ����������� ����������� ������� �� ������ � ��� (��������� ����������)
%> @warning ������������ ������� mapping �� ������� �����
% =========================================================================
    UpSempFlag = true(1);
    bits = randi([0 1], 1, 1000); % ��������� ���
    [sign,~] = mapping (bits, 'QPSK');       %QPSK 500 �������� 
    filtsign = filtration(sign, sqimpuls, nsamp, UpSempFlag);
    % =====================================================================
    %> �������� 3.1
    %> ��������� ������������ ������ �-�� � ������������������ �� ������������ ��������.
    % =====================================================================
    check3 = txfilter1(sign.').';
    if sum(abs(check3-filtsign))< 0.1 % �������� ���������� ���� 
                                      % ���������� ������������� 
                                      % � �������� ���������
        ans = '�������� ������ 3.1 �������� �������'
    else 
        err = '������ � ������ 3.1. ��������� ������'
        ans = sum(abs(check3-filtsign))
    end

    % ������ ����������� ���������
    figure;
    scatter(real(filtsign), imag(filtsign), 'o');
    title('���������� ��������� � ������������������');
    xlabel('In-Phase (I)');
    ylabel('Quadrature (Q)');
    grid on;
    
    % ��������� ������� �� �������
    figure;
    time = (0:length(filtsign)-1);
    plot(time, abs(filtsign));
    title('��������� ����������������������� �������');
    xlabel('�����');
    ylabel('���������');
    grid on;

    % =====================================================================
    %> �������� 3.2
    %> ��������� ������������ ������ �-�� ��� ����������������� �� ������������ ��������.
    % =====================================================================
    UpSempFlag = false(1);
    filtsign2 = filtration(filtsign, sqimpuls, nsamp, UpSempFlag);
    rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor', rolloff, ...
                                              'FilterSpanInSymbols',span,...
                                              'InputSamplesPerSymbol', nsamp,...
                                              'DecimationFactor', 1);
    chack4 = rxfilter(filtsign.').';
    if sum(abs(chack4-filtsign2))< 0.1 % �������� ���������� ���� 
                                       % ���������� ������������� 
                                       % � �������� ���������
        ans = '�������� ������ 3.2 �������� �������'
    else 
        err = '������ � ������ 3.2. ��������� ������'
        ans = sum(abs(chack4-filtsign2))
    end

     % ������ ����������� ���������
    figure;
    scatter(real(filtsign2), imag(filtsign2), 'o');
    title('���������� ��������� ��� �����������������');
    xlabel('In-Phase (I)');
    ylabel('Quadrature (Q)');
    grid on;
    
%     signal2 = filtsign2(2:nsamp:end);
%     scatterplot(signal2);
%     plot(abs(signal2));
    % ��������� ������� �� �������
    figure;
    time = (0:length(filtsign2)-1);
    plot(time, abs(filtsign2));
    title('��������� ������� ��� �����������������');
    xlabel('�����');
    ylabel('���������');
    grid on;


    %% point 5.1: screenshot model
    UpSempFlag = true(1);
    message = bits;
    [IQ,~] = mapping(message, 'QPSK');       %QPSK 500 �������� 
    IQ_filt = filtration(IQ, sqimpuls, nsamp, UpSempFlag);

    IQ_filt_noise=Noise(5,IQ_filt);

    IQ_filt2 = filtration(IQ_filt_noise, sqimpuls, nsamp, UpSempFlag);
    
    IQ_downsampled=downsample(IQ_filt2, nsamp);

    figure;
    scatter(real(IQ_filt2), imag(IQ_filt2), 'o');
    title('���������� ��������� �� �������� ����� ���������� AWGN');
    xlabel('In-Phase (I)');
    ylabel('Quadrature (Q)');
    grid on;
    
    figure;
    time = (0:length(IQ_filt2)-1);
    plot(time, abs(IQ_filt2));
    title('��������� ������� �� �������� ����� ���������� AWGN');
    xlabel('�����');
    ylabel('���������');
    grid on;

    figure; 
    scatter(real(IQ_downsampled), imag(IQ_downsampled), 'o');
    title('���������� ��������� �� �������� ����� ����������������� � �������� ������� �������������');
    xlabel('In-Phase (I)');
    ylabel('Quadrature (Q)');
    grid on;
    
    figure;
    time = (0:length(IQ_downsampled)-1);
    plot(time, abs(IQ_downsampled));
    title('��������� ������� �� �������� ����� ����������������� � �������� ������� �������������');
    xlabel('�����');
    ylabel('���������');
    grid on;
    
    %% point 5.2 MER(freq_offset)
    SNR=30; % in dB
    freqOffsetPercentage = -0.5:1/16:0.5;
    MER_values = zeros(size(freqOffsetPercentage));

    % ���� �� ������ ��������� ���������� ������
    for i = 1:length(freqOffsetPercentage)
        
        % �������� ��� � �������� SNR
        IQ_filt_noise=Noise(SNR,IQ_filt);

        % �������� ��������� �����
        frequencyOffset = freqOffsetPercentage(i);
        RX_IQ_offset = IQ_filt_noise .* exp(1i * 2 * pi * frequencyOffset * (0:length(IQ_filt_noise)-1));
        
        % ������������� ���������� �� ��������
        IQ_filt3 = filtration(RX_IQ_offset, sqimpuls, nsamp, 0);
        
        % ��������� �������������
        IQ_filt3_sync = IQ_filt3 .* exp(-1i * 2 * pi * frequencyOffset * (0:length(IQ_filt_noise)-1));

        % Downsampling �� ��������
        IQ_filt3_downsampled=downsample(IQ_filt3_sync, nsamp);

        % MER ��� �������� ���������� ������ 
        MER_values(i) = MER_my_func(IQ_filt3_downsampled, 'QPSK');
        %scatterplot (IQ_filt3_downsampled);
    end

    figure;
    plot(freqOffsetPercentage, MER_values, '-', 'LineWidth', 2);
    title('MER Characteristic vs. Frequency Offset');
    xlabel('Frequency Offset');
    ylabel('MER (dB)');
    grid on;