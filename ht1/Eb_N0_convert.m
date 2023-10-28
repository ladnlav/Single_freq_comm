function [Eb_N0] = Eb_N0_convert(SNR, Constellation)
    [~, bit_depth] = constellation_func(Constellation);
    Eb_N0 = SNR-10*log10(bit_depth);
end

