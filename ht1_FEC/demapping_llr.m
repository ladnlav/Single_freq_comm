function [LLR] = demapping_llr(pad, IQ, Constellation,N2_var,method)
% Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
% calculate the Bit_depth for each contellation
[Dictionary, Bit_depth_Dict] = constellation_func(Constellation);

% Amount of Const points
M=2^Bit_depth_Dict;

mapping=get_mapping(Constellation);
S_x = zeros(Bit_depth_Dict,M);

% Bit mapping matrix
for idx = 1:Bit_depth_Dict
    S_x(idx,:) =  bitget(mapping,idx);
end
S_x = logical(S_x);

% Find distances between a given IQ point and all points of Const
dist = zeros(2^Bit_depth_Dict,length(IQ));
for idx = 1:(2^Bit_depth_Dict)
    dist(idx,:) = (real(IQ) - (real(Dictionary(idx)))).^2 + (imag(IQ) - (imag(Dictionary(idx)))).^2;
end

if method == "exact"
    exp_mat = exp(-1./N2_var.*dist);
    % LLR calculation
    llr = zeros(Bit_depth_Dict,length(IQ));
    for idx = 1:Bit_depth_Dict
        llr(Bit_depth_Dict-idx+1,:) = log(sum(exp_mat(~S_x(idx,:),:),1)) - log(sum(exp_mat(S_x(idx,:),:),1));
    end
    LLR=llr(:);
elseif method == "app"
    llr_soft=zeros(Bit_depth_Dict,length(IQ));
    for idx = 1:Bit_depth_Dict
        llr_soft(Bit_depth_Dict-idx+1,:) = min(dist(~S_x(idx,:),:),[],1) - min(dist(S_x(idx,:),:),[],1);
    end
    llr_soft=llr_soft.*-1./N2_var;

    LLR=llr_soft(:);
else 
    disp("Unknown method!")
end

if pad~=-1
    LLR=LLR(1:end-pad);
end

end

function MAPPING = get_mapping(CONSTELLATION)
    % Define the mapping for different constellations
    switch CONSTELLATION
        case 'BPSK'
            MAPPING = [0; 1]; % BPSK mapping
        case 'QPSK'
            MAPPING = [0 1 2 3]; % QPSK mapping
        case '8PSK'
            MAPPING = [0 1 2 3 4 5 6 7]; % 8PSK mapping
        case '16QAM'
            MAPPING = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; % 16QAM mapping
        otherwise
            error('Unsupported constellation.');
    end
end
