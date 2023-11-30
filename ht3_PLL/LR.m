function [RX_IQ, LR_estimated] = LR(Channel_IQ)
    
LR_estimated = zeros(size(Channel_IQ(:,1)'));
RX_IQ = zeros(size(Channel_IQ));

for iter_time = 1: size(Channel_IQ,1)
    SOF = [1 0 0 1 1 1 0 1 0 1 0 1 0 1 1 0 0 1 0 0]; 
    IQ_SOF = mapping(SOF, 'BPSK'); % Use this sequence on the Rx as a Pilot-Signal
    N=19;

    % TASK
    r=Channel_IQ(iter_time,1:length(IQ_SOF));
    z=r.*conj(IQ_SOF);
    
    sum_angle=0;
    for m=1:N
        sum_angle=sum_angle+autocorr(m,z);
    end

    % LR detector
    LR_estimate = 1/(pi*(N+1))*angle(sum_angle);
    LR_estimated(iter_time)=LR_estimate;

    % TASK
    % NCO | Phase Accumulation
    LR_NCO = (1:size(Channel_IQ,2))*LR_estimate;
    % TASK    
    % Compensation
    RX_IQ(iter_time,:) = Channel_IQ(iter_time,:) .* exp(-1j * 2 * pi * LR_NCO);
end

% =========================================================================
% How does the estimate behave? Show on the plot
% How did the constellation change?
% -------------------------------------------------------------------------
end

function R=autocorr(m,z)
    
    L=length(z);
%     R=0;
%     for k=m+1:L
%         R=R+z(k)*conj(z(k-m));
%     end
% 
%     R=R/(L-m);

    R = 1 / (L - m) * sum(z(m+1:L) .* conj(z(1:L-m)));
end