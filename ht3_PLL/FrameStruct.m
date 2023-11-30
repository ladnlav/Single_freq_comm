function [Data_By_Frame_With_Header] = FrameStruct(Data, Header, Amount_of_Frame)
    
    Data_By_Frame = reshape(Data, Amount_of_Frame, []);
    Copy_header = repmat(Header,Amount_of_Frame, 1);
    Data_By_Frame_With_Header = [Copy_header, Data_By_Frame];

   % Data_By_Frame_With_Header=reshape(Data_By_Frame_With_Header.',1,[]);

end

