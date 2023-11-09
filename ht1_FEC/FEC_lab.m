clc; clear;
%% ============================= Part 1  Warmup =====================================
% configuring LDPC encoder and decoder

% prototype matrix as defnied in Wi-Fi (IEEE® 802.11)
P = [
    16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
    25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
    25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
     9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
    24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
     2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 27; % N.B. this blocksize is connected with optimized pairty check matrix generation method 
                % it is NOT the blocksize of the ldpc
                % https://prezi.com/aqckvai6jux-/ldpc/?utm_campaign=share&utm_medium=copy
H = ldpcQuasiCyclicMatrix(blockSize,P); % getting parity-check matrix

cfgLDPCEnc = ldpcEncoderConfig(H); % configuring encoder
cfgLDPCDec = ldpcDecoderConfig(H); % configuring decoder

% using cfgLDPCEnc variable, print our the number of inofrmation, parity
% check bits and the coderate
fprintf('Number of information bits in a block: %d\n', cfgLDPCEnc.NumInformationBits);
fprintf('Number of parity check bits in a block: %d\n', cfgLDPCEnc.NumParityCheckBits);
coderate = cfgLDPCEnc.NumInformationBits / cfgLDPCEnc.BlockLength;
fprintf('Coderate: %f\n', coderate);

% simple test to check that encoder and decoder configured correctly
test_message = boolean(randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8'));
encodedData = ldpcEncode(test_message,cfgLDPCEnc);

% calculate the syndrome
s = encodedData'*H';
s=mod(s, 2); % we need xor instead of multiplication
if(~any(s))
    fprintf('No errors!\n');
else
    fprintf('Errors detected during syndrome check!\n');
end

% deliberately distorting one bit of the message
flip_index=randi(numel(encodedData));
encodedData(flip_index) = ~(encodedData(flip_index));

% checking the syndrome once again
s = encodedData'*H';
s=mod(s, 2); % we need xor instead of multiplication
if(~any(s))
    fprintf('No errors!\n');
else
    fprintf('Errors detected during syndrome check!\n');
end

%% ============= Part 2 comparing coded and uncoded system =================
Length_Bit_vector = 1e6;
rng(321); % Fix the seed of the random number generator
Constellation = "16QAM"; % BPSK, QPSK, 8PSK, 16QAM

% % Bit generator
% Bit_Tx = generateBits(Constellation,Length_Bit_vector);
% % Mapping
% IQ_TX = mapping(Bit_Tx, Constellation);
% % Eb_N0_convert(), which convert SNR to Eb/N0
% Eb_N0 = Eb_N0_convert(SNR, Constellation);
% % adding AWGN
% IQ_RX = Noise(SNR, IQ_TX);
% scatterplot(IQ_RX);
% % Demapping
% Bit_Rx = demapping(IQ_RX, Constellation);
% % Error check
% % Write your own function Error_check() for calculation of BER
% BER = Error_check(Bit_Tx, Bit_Rx);


%%
maxnumiter = 10;
snr = 0:22;
numframes = 10000;

% check manual on the build-in ber counter
% it outputs three variables
ber = comm.ErrorRate; %build-in BER counter
ber2 = comm.ErrorRate; %build-in BER counter
ber3 = comm.ErrorRate; %build-in BER counter

% arrays to store error statistic
errStats = zeros(length(snr), numframes, 3);
errStats_app = zeros(length(snr), numframes, 3); 
errStatsNoCoding = zeros(length(snr), numframes, 3);

% times for app and axact llr
t_exact = zeros(length(snr), 1);
t_app = zeros(length(snr), 1);

tStart = tic;
for ii = 1:length(snr)
    t_exact_i = 0;
    t_app_i = 0;
    parfor counter = 1:numframes
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
        % Transmit and receive with LDPC coding
        encodedData = ldpcEncode(data,cfgLDPCEnc);
        
        % YOUR MAPPER HERE choose any constellation type you like
        [modSignal,pad] = mapping(encodedData, Constellation);

        [rxsig, noisevar] = Noise(snr(ii),modSignal); % use yiur AWGN function

        % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
        t_exact_start = tic;
        llr_exact = demapping_llr(pad,rxsig, Constellation, noisevar^2,'exact');
        t_exact_i = t_exact_i + toc(t_exact_start);

        t_app_start = tic;
        llr_app = demapping_llr(pad,rxsig, Constellation, noisevar^2,'app');
        t_app_i = t_app_i + toc(t_app_start);

        rxbits = ldpcDecode(llr_exact,cfgLDPCDec,maxnumiter);
        errStats(ii, counter, :) = ber(data,rxbits);

        rxbits = ldpcDecode(llr_app,cfgLDPCDec,maxnumiter);
        errStats_app(ii, counter, :) = ber3(data,rxbits);
        %========================================
        
        % no coding system
        [noCoding, pad] = mapping(data, Constellation);
        [rxNoCoding, noisevar] = Noise(snr(ii),noCoding); % use your AWGN function
        % YOUR DEMAPPER HERE N.B. Perform Hard Demapping, output bits!
        rxBitsNoCoding = demapping(pad, rxNoCoding, Constellation);
        errStatsNoCoding(ii, counter, :) = ber2(data,int8(rxBitsNoCoding'));
    end

    % Calculate average time for this iteration
    t_exact(ii) = t_exact_i / numframes;
    t_app(ii) = t_app_i / numframes;

    fprintf(['SNR = %2d\n   Coded: Error rate = %1.5f, ' ...
        'Number of errors = %d\n'], ...
        snr(ii),mean(errStats(ii, :, 1), 2), mean(errStats(ii, :, 2), 2))
    fprintf(['Noncoded: Error rate = %1.5f, ' ...
        'Number of errors = %d\n'], ...
        mean(errStatsNoCoding(ii, :, 1), 2), mean(errStatsNoCoding(ii, :, 2), 2))
    reset(ber);
    reset(ber2);
    reset(ber3);
end

ber.release();
ber2.release();
tend = toc(tStart);
fprintf('Simulation finished after %.2f s\n', tend);
%%
%
figure();
semilogy(snr, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2)
hold on
semilogy(snr, mean(errStats(:, :, 1), 2), 'LineWidth', 2)
hold off

xlabel('SNR, dB');
ylabel('BER')
grid on
set(gca, 'Fontsize', 20)
legend("NoCoding", "LDPC",'Location','southwest');
title('Зависимость BER от SNR с LDPC и без');

save('BER_SNR_results.mat', 'errStatsNoCoding', 'errStats', '-v7.3')
% 
% Replot the results in BER vs Eb/N0 scale

Eb_N0=Eb_N0_convert(snr,Constellation);
figure();
semilogy(Eb_N0, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2)
hold on
semilogy(Eb_N0, mean(errStats(:, :, 1), 2), 'LineWidth', 2)
hold off

xlabel('Eb/N0, dB');
ylabel('BER')
grid on
set(gca, 'Fontsize', 20)
legend("NoCoding", "LDPC",'Location','southwest')
title('Зависимость BER от Eb/N0 с LDPC и без')
% how the shape of curves has changed?

% The shape is the same. The curves moved to the left along the x axis and
% that's it.

% % what is the gain in dB?
M_nocod=mean(errStatsNoCoding(:, :, 1), 2);
M_code=mean(errStats(:, :, 1), 2);
Gain=(10*log10((M_nocod(M_code~=0)./M_code(M_code~=0))));
Gain_max=max(Gain);
disp(Gain_max);

figure();
plot(Eb_N0(M_code~=0), Gain, 'LineWidth', 2);

xlabel('Eb/N0, dB');
ylabel('GAIN, dB')
grid on
set(gca, 'Fontsize', 20)
title('Gain в зависимости от SNR')

% +20 points: compare results with llr and approximate llr formulas
figure();
semilogy(Eb_N0, mean(errStats(:, :, 1), 2), 'LineWidth', 2)
hold on
semilogy(Eb_N0, mean(errStats_app(:, :, 1), 2),'--', 'LineWidth', 2)
hold off

xlabel('Eb/N0, dB');
ylabel('BER')
grid on
set(gca, 'Fontsize', 20)
legend("Exact", "Approximate",'Location','southwest')
title('Сравнение результатов работы точного и приблизительного расчетов LLR')

% Time comparison
figure();
plot(Eb_N0, t_exact', 'LineWidth', 2);
hold on;
plot(Eb_N0,t_app', '--', 'LineWidth', 2);
xlabel('Eb/N0, dB');
ylabel('Time, sec')
grid on
set(gca, 'Fontsize', 20)
legend('Exact', 'Approximate','Location','southwest')
title('Время работы двух вариантов расчёта LLR')

%% ================ Part 3: default LDPC with different numbers of iterations =========================

%change the snr range to capture behaviour of coded curves only
snr2 = 6:0.2:14;

maxnumiters = [5, 10, 20, 30]; % we will plot curves for two values of decoding iterations
numframes = 5000;
errStats_it_num = zeros(length(snr2), numframes, 3, numel(maxnumiters));

tStart = tic;

% +10 points for using parfor here and calculating speedup
parfor ii = 1:length(snr2)
    stats=zeros(numframes, 3, numel(maxnumiters));
    for m = 1:numel(maxnumiters)
        maxnumiter = maxnumiters(m);
        for counter = 1:numframes
            data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
            % Transmit and receive with LDPC coding
            encodedData = ldpcEncode(data,cfgLDPCEnc);

            % YOUR MAPPER HERE choose any constellation type you like
            [modSignal,pad] = mapping(encodedData,Constellation);

            [rxsig, noisevar] = Noise(snr2(ii),modSignal); % use yiur AWGN function

            % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
            llr = demapping_llr(pad, rxsig,Constellation,noisevar^2,'app');

            rxbits = ldpcDecode(llr,cfgLDPCDec,maxnumiter);
            br = ber(data,rxbits);
            stats(counter, :, m) = br;
        end
%         fprintf(['SNR = %2d\n   Coded with %d iterations: Error rate = %1.5f, ' ...
%             'Number of errors = %d\n'], ...
%             snr2(ii), maxnumiter, mean(errStats_it_num(ii, :, 1, m), 2), mean(errStats_it_num(ii, :, 2, m), 2))
        reset(ber);
    end
    errStats_it_num(ii,:,:,:)=stats;
end
ber.release();
tend = toc(tStart);
fprintf('Simulation finished in %.2f s\n', tend);

%%
Eb_N0=Eb_N0_convert(snr,Constellation);
Eb_N02=Eb_N0_convert(snr2,Constellation);
figure();
semilogy(Eb_N0, mean(errStatsNoCoding(:, :, 1), 2), 'LineWidth', 2)
hold on
for m = 1:numel(maxnumiters)
    semilogy(Eb_N02, mean(errStats_it_num(:, :, 1, m), 2), 'LineWidth', 2)
    hold on
end
hold off
grid on
xlabel('Eb/N0, dB')
ylabel('BER')
set(gca, 'FontSize', 20)

legend('No coding', strcat('LDPC 3/4 ', 32, num2str(maxnumiters(1)), 32,'iterations'), ...
    strcat('LDPC 3/4 ', 32, num2str(maxnumiters(2)), 32, 'iterations'), ...
    strcat('LDPC 3/4 ', 32, num2str(maxnumiters(3)), 32,'iterations'), ...
    strcat('LDPC 3/4 ', 32, num2str(maxnumiters(4)), 32,'iterations'),'Location','southwest');
% 
save('BER_SNR_results.mat', 'errStats_it_num', '-append')

% change the plot to Eb/N0 scale! - DONE
%% ========================= Part 4: diffrent decoding methods with the same max number of iterations

% https://www.mathworks.com/help/comm/ref/ldpcdecode.html
cfgLDPCDec2 = ldpcDecoderConfig(H, 'norm-min-sum'); % configuring second decoder
maxnumiter = 10;
snr2 = 6:0.2:14;
numframes = 10000;


errStats_minsum = zeros(length(snr2), numframes, 3);
errStats_bp = zeros(length(snr2), numframes, 3);

ber = comm.ErrorRate;
ber2 = comm.ErrorRate;

MinSumScalingFactor_array=0.1:0.1:1;
%MinSumScalingFactor = 0.75; % task: find the best parameter 
t_min_sum = zeros(length(snr2), 1);
t_bp = zeros(length(snr2), 1);
errStats_minsum_array=zeros(numel(MinSumScalingFactor_array),length(snr2), numframes, 3);

tStart = tic;

for jj=1:numel(MinSumScalingFactor_array)
MinSumScalingFactor = MinSumScalingFactor_array(jj);
for ii = 1:length(snr2)
    t_min_sum_i = 0;
    t_bp_i = 0;
    parfor counter = 1:numframes
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
        % Decode with belief propagation
        encodedData = ldpcEncode(data,cfgLDPCEnc);

        % YOUR MAPPER HERE choose any constellation type you like
        [modSignal,pad] = mapping(encodedData,Constellation);

        [rxsig, noisevar] = Noise(snr2(ii),modSignal); % use yiur AWGN function
        
        % YOUR DEMAPPER HERE N.B. Perform Soft Demapping, output llr!
        llr = demapping_llr(pad, rxsig,Constellation,noisevar^2,'app');
        
        % decode with MinSum
        t_min_sum_start = tic;
        rxbits = ldpcDecode(llr,cfgLDPCDec2,maxnumiter, 'MinSumScalingFactor', MinSumScalingFactor);
        t_min_sum_i = t_min_sum_i + toc(t_min_sum_start);

        errStats_minsum(ii, counter, :) = ber(data,rxbits);

        % ================================
        % Decode with layered belief propagation
        t_bp_start = tic;
        rxbits = ldpcDecode(llr,cfgLDPCDec,maxnumiter);
        t_bp_i = t_bp_i + toc(t_bp_start);

        errStats_bp(ii, counter, :) = ber2(data,rxbits);
    end
    
    % Calculate average time for this iteration
    t_min_sum(ii) = t_min_sum_i / numframes;
    t_bp(ii) = t_bp_i / numframes;

    fprintf(['SNR = %2d\n   Min Sum decoding: Error rate = %e, ' ...
        'Number of errors = %d, average time %.4f s\n'], ...
        snr2(ii),mean(errStats_minsum(ii, :, 1), 2), mean(errStats_minsum(ii, :, 2), 2),  t_min_sum(ii))
    fprintf(['BP decoding: Error rate = %e, ' ...
        'Number of errors = %d, average time %.4f s\n'], ...
        mean(errStats_bp(ii, :, 1), 2), mean(errStats_bp(ii, :, 2), 2), t_bp(ii))
    reset(ber);
    reset(ber2);
end

errStats_minsum_array(jj,:,:,:)=errStats_minsum;
end
t = toc(tStart);
fprintf('Simulation finished after %.2f s\n', t);
%% MinSum algorithm with different MinSumScalingFactor
figure();
Eb_N02=Eb_N0_convert(snr2,Constellation);
legendLabels = cell(1, numel(MinSumScalingFactor_array));

for i=1:numel(MinSumScalingFactor_array)
    
    semilogy(Eb_N02, mean(errStats_minsum_array(i,:, :, 1), 3),'LineWidth', 2)
    hold on
    legendLabels{i}=['MinSumScalingFactor = ' num2str(MinSumScalingFactor_array(i))];
end

legend(legendLabels,'Location','southwest')
xlabel('E_b/N_0, dB')
ylabel('BER')
grid on
set(gca, 'FontSize', 20)
title('MinSum алгоритм с разными MinSumScalingFactor')
%%
Eb_N02=Eb_N0_convert(snr2,Constellation);
semilogy(Eb_N02, mean(errStats_minsum(:, :, 1), 2))
hold on
semilogy(Eb_N02, mean(errStats_bp(:, :, 1), 2), '--')
hold off

legend('MinSum', 'Belief Propagation','Location','southwest')
xlabel('E_b/N_0, dB')
ylabel('BER')
grid on
set(gca, 'FontSize', 20)

save('BER_SNR_results.mat', 'errStats_minsum', 'errStats_bp', '-append')

%% Part four: compare the speed of the algorithms
% compare the speed of Belief Propagation and MinSum decoders

figure();
plot(Eb_N02, t_min_sum', 'LineWidth', 2);
hold on;
plot(Eb_N02,t_bp', '--', 'LineWidth', 2);
xlabel('E_b/N_0, dB');
ylabel('Time, sec')
grid on
set(gca, 'Fontsize', 20)
legend('MinSum', 'Belief Propagation','Location','southwest')
title('Время работы двух алгоритмов LDPC')