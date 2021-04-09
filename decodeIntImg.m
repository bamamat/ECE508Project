function [RGBrecovered, metaData] = decodeIntImg(binaryReceived, txAntCalc, nTxAnt)
    txAntCalcBin = dec2bin(txAntCalc, log2(nTxAnt));
    data = num2str(binaryReceived);
    data(isspace(data)) = '';
    data = [txAntCalcBin data];

    metaLen = bin2dec(data(1:8));
    meta = zeros(1,metaLen);
    for metaCount = 1:metaLen
        idx = 9 + 16*(metaCount-1);
        meta(metaCount) = bin2dec(data(idx:idx+15));
    end

    data = data(idx+16:end);
    RGBout = uint8(zeros(1,length(data)/8));
    for dataCount = 1:length(data)/8
        idx = 1 + 8*(dataCount-1);
        RGBout(dataCount) = uint8(bin2dec(data(idx:idx+7)));
    end
    RGBrecovered = reshape(RGBout, meta(1),meta(2),meta(3));
    metaData = meta(4:end);
end