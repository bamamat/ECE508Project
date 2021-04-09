function [binaryIntVect] = encodeIntImg(RGBin, int16MetaData)
    metaData = [size(RGBin) int16MetaData];
    RGB = reshape(RGBin,1,numel(RGBin));
    dataToSend = [dec2bin(length(metaData),8) reshape(dec2bin(metaData,16).',1,[]) reshape(dec2bin(uint8(RGB),8).',1,[])];
    binaryIntVect = zeros(1,length(dataToSend));
    for dataIdx = 1:length(dataToSend)
        if dataToSend(dataIdx) == '1'
            binaryIntVect(dataIdx) = 1;
        else
            binaryIntVect(dataIdx) = 0;
        end
    end
end