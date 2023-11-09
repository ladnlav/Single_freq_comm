function bits = generateBits(constellation,Length_Bit_vector)

    [~, Bit_depth_Dict] = constellation_func(constellation);

    bits = randi([0 1], 1, Length_Bit_vector);

    % Calculate the number of zeros needed to pad the bits vector
    num_zeros = Bit_depth_Dict - mod(length(bits), Bit_depth_Dict);
    if num_zeros == Bit_depth_Dict
        num_zeros = 0;
    end
    
    % Pad the bits vector with zeros
    bits = [bits zeros(1, num_zeros)];
end