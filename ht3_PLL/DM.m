function [RX_IQ_DM, DM_estimate, DM_Filtred_ar] = DM(Channel_IQ, Kp, Ki)
 % Initialization
    DM_NCO = 0;  % Initial phase of the NCO
    DM_estimate = zeros(size(Channel_IQ(:,1)'));  % Initialize DM estimate
    RX_IQ_DM = zeros(size(Channel_IQ));
    DM_Filtred_ar=zeros(size(DM_estimate));

    integrator = zeros(1,numel(DM_estimate)+1);
 % Parameters
    D = 2;  % Delay factor
    SOF = [1 0 0 1 1 1 0 1 0 1 0 1 0 1 1 0 0 1 0 0]; 
    IQ_SOF = mapping(SOF, 'BPSK'); % Use this sequence on the Rx as a Pilot-Signal

   
 % Loop processing
    for iter_time = 1:size(Channel_IQ,1)
        % Compensation
        Channel_IQ(iter_time,:)=Channel_IQ(iter_time,:) .* exp(-1j * 2 * pi * DM_NCO);
        RX_IQ = Channel_IQ(iter_time,:);

        % DM detector
        r=Channel_IQ(iter_time,1:length(IQ_SOF));
    
        z=r.*conj(IQ_SOF);
%         sum_df=0;
%         for k = D+1:length(IQ_SOF)
%             sum_df=sum_df+z(k)*conj(z(k-D));
%         end
        sum_df=sum(z(D+1:end).*conj(z(1:end-D)));
        DM_estimate(iter_time) = 1/(2*pi*D)*angle(sum_df);
        
        % Loop filter
        integrator(iter_time+1)=Ki*DM_estimate(iter_time)+integrator(iter_time);
        DM_Filtred = integrator(iter_time+1) + Kp*DM_estimate(iter_time);
        

        % NCO | Phase Accumulation
        DM_NCO = (1:size(Channel_IQ,2))*DM_Filtred;

        % Output
        RX_IQ_DM(iter_time,:) = RX_IQ;
        DM_Filtred_ar(iter_time)=DM_Filtred;
    end
    
   %scatter(real(RX_IQ_DM),imag(RX_IQ_DM))
   
end
% =========================================================================
% TASK
% For different Damping Factor and BnTs calculate coefficients of loop filter
% What changes in synchronisation when the loop filter coefficients are recalculated?
% Illustrate these changes on the graphs
% How did the constellation change?
% -------------------------------------------------------------------------


